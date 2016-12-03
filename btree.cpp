/**
 * @author See Contributors.txt for code contributors and overview of BadgerDB.
 *
 * @section LICENSE
 * Copyright (c) 2012 Database Group, Computer Sciences Department, University of Wisconsin-Madison.
 */

#include "btree.h"
#include "filescan.h"
#include <assert.h>
#include "exceptions/bad_index_info_exception.h"
#include "exceptions/bad_opcodes_exception.h"
#include "exceptions/bad_scanrange_exception.h"
#include "exceptions/no_such_key_found_exception.h"
#include "exceptions/scan_not_initialized_exception.h"
#include "exceptions/index_scan_completed_exception.h"
#include "exceptions/file_not_found_exception.h"
#include "exceptions/end_of_file_exception.h"

#define EMPTY_SLOT ((1)<<(sizeof(int)*8 - 1))

//#define DEBUG

namespace badgerdb
{

	// -----------------------------------------------------------------------------
	// BTreeIndex::BTreeIndex -- Constructor
	// -----------------------------------------------------------------------------

	BTreeIndex::BTreeIndex(const std::string & relationName,
			std::string & outIndexName,
			BufMgr *bufMgrIn,
			const int attrByteOffset,
			const Datatype attrType)
	{
		// Constructing an index name 
		std::ostringstream idxStr;
		idxStr << relationName << '.' << attrByteOffset;
		outIndexName = idxStr.str(); // indeName is the name of the index file
		std::cout << "empty slot is" << EMPTY_SLOT << std::endl; // output test for the the opened file
		std::cout << "Name of index is: " << outIndexName << std::endl; // output test for the the opened file

		// Declare a page instance
		Page *page;

		// Define the value of bufMgr
		bufMgr = bufMgrIn;

		// If the index file exists, the file is opened
		// Check whether the file opened matches the buffer info
		try
		{

			file = new BlobFile(outIndexName, false);
			assert(file != NULL);
			std::cout << "The file is opened." << std::endl; // output test for the the opened file

			// Read the header page of from the buffer pool
			bufMgr -> readPage (file, 1, page); 

			// Cast page to IndexMetaInfo to retrieve info
			IndexMetaInfo* indexMetaInfo = reinterpret_cast<IndexMetaInfo*>(page);

			// Unpin header page since it is not required
			bufMgr -> unPinPage(file, 1, false);

			if (!_validateMetaPage(indexMetaInfo, relationName, attrByteOffset, attrType)) {
				throw BadIndexInfoException(outIndexName);
			}

			
			rootPageNum = indexMetaInfo->rootPageNo;
		}

		// If the index file does not exist, then a new file is created
		catch (FileNotFoundException e)
		{
			// If the file does not exist, a new index file is created
			file = new BlobFile(outIndexName, true);
			std::cout << "A new index file is created";// Test output message    

			{
				PageId pid;
				IndexMetaInfo *meta;
				bufMgr->allocPage(file, pid, (Page *&)meta);
				strncpy((char *)(&(meta->relationName)), relationName.c_str(), 20);
				meta->relationName[19] = 0;
				meta->attrByteOffset = attrByteOffset;
				meta->attrType = attrType;

				{
					PageIDPair rootResult = _newNonLeafNode();
					rootPageNum = rootResult.pageNo;
					NonLeafNodeInt *root = (NonLeafNodeInt *)(rootResult.page);
					root->level = 1;

					PageIDPair leafResult = _newLeafNode();
					bufMgr->unPinPage(file, leafResult.pageNo, true);
					root->pageNoArray[0] = leafResult.pageNo;

					bufMgr->unPinPage(file, rootResult.pageNo, true);
				}

				meta->rootPageNo = rootPageNum;
				bufMgr->unPinPage(file, pid, true);
			}
			

			// Scan the relationship(using FileScan) 
			FileScan *myFileScan = new FileScan(relationName, bufMgrIn);

			// Scan the relationship using FileScan and insert the entries for all the tuples
			// in this relation into the index
			try 
			{
				while (1)
				{
					RecordId outRid;
					myFileScan -> scanNext(outRid);

					// Using getRecord() method to get all the record in the file
					std::string record = myFileScan -> getRecord();
					const char *cstr = record.c_str();
					// std::cout << "My record: " << record << std::endl;
					// std::cout << outRid.page_number << "   " << outRid.slot_number << "\n";
					insertEntry(cstr + attrByteOffset, outRid);
				}
			}
			catch (EndOfFileException e)
			{
				std::cout << "Reach the end of the file." << std::endl;
			}
			delete myFileScan;
			
		}


		// Test for testing the file is opened
		assert(file->isOpen(outIndexName));    


	}

	const bool BTreeIndex::_validateMetaPage(IndexMetaInfo *meta, std::string relationName, int attrByteOffset, int attrType)
	{
		// If not match, throw bad_index_info_exception
		if (relationName != meta -> relationName || attrByteOffset != meta -> 
				attrByteOffset || attrType != meta -> attrType)
		{
			return false;
		}
		return true;
	}

	const RIDKeyPair<int> BTreeIndex::_getRIDKeyPairFromRecord(std::string record, RecordId rid, int offset)
	{
		RIDKeyPair<int> pair;
		const char *cstr = record.c_str();
		int key;
		key = *(cstr+offset);
		pair.set(rid, key);
		return pair;	
	}

	// -----------------------------------------------------------------------------
	// BTreeIndex::~BTreeIndex -- destructor
	// -----------------------------------------------------------------------------

	BTreeIndex::~BTreeIndex()
	{
		bufMgr->flushFile(file);
                //delete file;
		file->~File();
	}

	// -----------------------------------------------------------------------------
	// BTreeIndex::insertEntry
	// -----------------------------------------------------------------------------

	const void BTreeIndex::insertEntry(const void *key, const RecordId rid) 
	{
		PageKeyPair<int> result = _insertEntryNonLeaf(rootPageNum, key, rid);
		if (result.pageNo != (PageId)(EMPTY_SLOT)) {
			PageIDPair newRootPair = _newNonLeafNode();
			NonLeafNodeInt *newRoot = (NonLeafNodeInt *)(newRootPair.page);
			newRoot->level = 0;
			newRoot->keyArray[0] = result.key;
			newRoot->pageNoArray[0] = rootPageNum;
			newRoot->pageNoArray[1] = result.pageNo;
			bufMgr->unPinPage(file, newRootPair.pageNo, true);

			rootPageNum = newRootPair.pageNo;

			IndexMetaInfo *meta;
			bufMgr->readPage(file, 1, (Page *&)meta);
			meta->rootPageNo = rootPageNum;
			bufMgr->unPinPage(file, 1, true);
		}
	}

	const PageKeyPair<int> BTreeIndex::_insertEntryNonLeaf(PageId nodeId, const void* key, const RecordId rid)
	{
		NonLeafNodeInt *node;
		bufMgr->readPage(file, nodeId, (Page *&)node);
		PageKeyPair<int> retVal;
		retVal.set(EMPTY_SLOT, EMPTY_SLOT);
		bool dirty = false;
		{
			int ikey = *((int *)key);
			int idx = _search(node, ikey);
			PageId nextNodeId = node->pageNoArray[idx];
			PageKeyPair<int> result;
			if (node->level == 1) {
				result = _insertEntryLeaf(nextNodeId, key, rid);
			} else {
				result = _insertEntryNonLeaf(nextNodeId, key, rid);
			}
			if (result.pageNo == PageId(EMPTY_SLOT)) {
				goto done;
			}
			if (!_nonLeafIsFull(node)) {
				_nonLeafInsertEntry(node, result);
			} else {
				retVal = _nonLeafSplitInsertEntry(node, result);
			}
			dirty = true;
		}
done:
		bufMgr->unPinPage(file, nodeId, dirty);
		return retVal;
	}

	const PageKeyPair<int> BTreeIndex::_insertEntryLeaf(PageId nodeId, const void* key, const RecordId rid)
	{
		LeafNodeInt *node;
		bufMgr->readPage(file, nodeId, (Page *&)node);
		PageKeyPair<int> retVal;
		retVal.set(EMPTY_SLOT, EMPTY_SLOT);
		bool dirty = false;
		{
			if (!_leafIsFull(node)) {
				_leafInsertEntry(node, key, rid);
			} else {
				retVal = _leafSplitInsertEntry(node, key, rid);
			}
			dirty = true;
		}
		bufMgr->unPinPage(file, nodeId, dirty);
		return retVal;
	}

	const bool BTreeIndex::_leafIsFull(LeafNodeInt *node)
	{
		return node->keyArray[INTARRAYLEAFSIZE - 1] != EMPTY_SLOT;
	}

	const bool BTreeIndex::_nonLeafIsFull(NonLeafNodeInt *node)
	{
		return node->keyArray[INTARRAYNONLEAFSIZE - 1] != EMPTY_SLOT;
	}

	const void BTreeIndex::_leafInsertEntry(LeafNodeInt *node, const void *key, const RecordId rid)
	{
		assert(!_leafIsFull(node));

		int ikey = *((int *)key);
		//assert(ikey >= node->keyArray[0]);

		int emptySlotIdx = -1;
		for (int i = 0; i < INTARRAYLEAFSIZE; i++) {
			if (node->keyArray[i] == EMPTY_SLOT) {
				emptySlotIdx = i;
				break;
			}
		}

		for (int i = emptySlotIdx; i >= 0; i--) {
			if (i == 0) {
				node->keyArray[0] = ikey;
				node->ridArray[0] = rid;
			} else {
				if (node->keyArray[i-1] > ikey) {
					node->keyArray[i] = node->keyArray[i-1];
					node->ridArray[i] = node->ridArray[i-1];
				} else {
					node->keyArray[i] = ikey;
					node->ridArray[i] = rid;
					break;
				}
			}
		}
		_assertLeafInternalConsistency(node);
	}

	const void BTreeIndex::_nonLeafInsertEntry(NonLeafNodeInt *node, const PageKeyPair<int> pair)
	{
		assert(!_nonLeafIsFull(node));

		int emptySlotIdx = -1;
		for (int i = 0; i < INTARRAYNONLEAFSIZE; i++) {
			if (node->keyArray[i] == EMPTY_SLOT) {
				emptySlotIdx = i;
				break;
			}
		}
		
		for (int i = emptySlotIdx; i >= 0; i--) {
			if (i == 0) {
				node->keyArray[0] = pair.key;
				node->pageNoArray[1] = pair.pageNo;
			} else {
				if (node->keyArray[i-1] > pair.key) {
					node->keyArray[i] = node->keyArray[i-1];
					node->pageNoArray[i+1] = node->pageNoArray[i];
				} else {
					node->keyArray[i] = pair.key;
					node->pageNoArray[i+1] = pair.pageNo;
					break;
				}
			}
		}
		_assertNonLeafInternalConsistency(node);
	}

	const PageKeyPair<int> BTreeIndex::_leafSplitInsertEntry(LeafNodeInt *node, const void *key,  RecordId rid) {
		assert(_leafIsFull(node));
		int half = (INTARRAYLEAFSIZE  + 1)/ 2;
		PageIDPair newLeaf =_newLeafNode();
		LeafNodeInt *newNode = reinterpret_cast<LeafNodeInt *>(newLeaf.page);

		int ikey = *((int *)key);

		bool insertLeft = false;
		if (ikey < node->keyArray[half])
		{
			insertLeft = true;
			half = half - 1;
		}

		for (int i = half; i < INTARRAYLEAFSIZE; i++)
		{
			newNode->keyArray[i - half] = node->keyArray[i];
			newNode->ridArray[i - half] = node->ridArray[i];
			node->keyArray[i] = EMPTY_SLOT;
			node->ridArray[i].page_number = (PageId) EMPTY_SLOT;
			node->ridArray[i].slot_number = (SlotId) EMPTY_SLOT;
		}
		if (insertLeft) 
		{
			_leafInsertEntry(node, key, rid);

		} else {	
			_leafInsertEntry(newNode, key, rid);
		}	
		PageKeyPair<int> newPair;
		newPair.set(newLeaf.pageNo, newNode->keyArray[0]); 
		newNode->rightSibPageNo = node->rightSibPageNo;
		node->rightSibPageNo = newLeaf.pageNo;

		_assertLeafInternalConsistency(node);
		_assertLeafInternalConsistency(newNode);

		bufMgr->unPinPage(file, newPair.pageNo, true);
		return newPair;
	}

	const PageKeyPair<int> BTreeIndex::_nonLeafSplitInsertEntry(NonLeafNodeInt *node, const PageKeyPair<int> pair)
	{
		PageIDPair newSibPair =_newNonLeafNode();
		NonLeafNodeInt *newSib = reinterpret_cast<NonLeafNodeInt *>(newSibPair.page);

		newSib->level = node->level;
		PageKeyPair<int> parentEntry;

		int index = 0;
		for (; index < INTARRAYNONLEAFSIZE && node->keyArray[index] < pair.key; index++);

		int half = (INTARRAYNONLEAFSIZE + 1) / 2;
		if (index < half) {
			for (int i = 0; i < half - 1; i++) {
				newSib->keyArray[i] = node->keyArray[i + half];
				newSib->pageNoArray[i] = node->pageNoArray[i + half];
				node->keyArray[i + half] = EMPTY_SLOT;
				node->pageNoArray[i + half] = EMPTY_SLOT;
			}
			newSib->pageNoArray[half - 1] = node->pageNoArray[INTARRAYNONLEAFSIZE];
			node->pageNoArray[INTARRAYNONLEAFSIZE] = EMPTY_SLOT;

			parentEntry.set(newSibPair.pageNo, node->keyArray[half - 1]);

			node->keyArray[half - 1] = EMPTY_SLOT;
			_nonLeafInsertEntry(node, pair);
		} else if (index == half) {
			for (int i = 0; i < half - 1; i++) {
				newSib->keyArray[i] = node->keyArray[i + half];
				newSib->pageNoArray[i+1] = node->pageNoArray[i + half + 1];
				node->keyArray[i + half] = EMPTY_SLOT;
				node->pageNoArray[i + half + 1] = EMPTY_SLOT;
			}
			newSib->pageNoArray[0] = pair.pageNo;

			parentEntry.set(newSibPair.pageNo, pair.key);
		} else {
			for (int i = 0; i < half - 2; i++) {
				newSib->keyArray[i] = node->keyArray[i + half + 1];
				newSib->pageNoArray[i] = node->pageNoArray[i + half + 1];
				node->keyArray[i + half + 1] = EMPTY_SLOT;
				node->pageNoArray[i + half + 1] = EMPTY_SLOT;
			}
			newSib->pageNoArray[half - 2] = node->pageNoArray[INTARRAYNONLEAFSIZE];
			node->pageNoArray[INTARRAYNONLEAFSIZE] = EMPTY_SLOT;

			parentEntry.set(newSibPair.pageNo, node->keyArray[half]);

			node->keyArray[half] = EMPTY_SLOT;
			_nonLeafInsertEntry(newSib, pair);
		}

		_assertNonLeafInternalConsistency(node);
		_assertNonLeafInternalConsistency(newSib);

		bufMgr->unPinPage(file, newSibPair.pageNo, true);
		return parentEntry;
	}

	// -----------------------------------------------------------------------------
	// BTreeIndex::startScan
	// -----------------------------------------------------------------------------

	const void BTreeIndex::startScan(const void* lowValParm,
			const Operator lowOpParm,
			const void* highValParm,
			const Operator highOpParm)
	{
		lowValInt = *((int *)lowValParm);
		highValInt = *((int *)highValParm);
		lowOp = lowOpParm;
		highOp = highOpParm;
		_checkScanParams(lowValInt, lowOpParm, highValInt, highOpParm);

		NonLeafNodeInt *node;
		LeafNodeInt *leafNode;
		PageId pageNum = rootPageNum;

		bufMgr->readPage(file, pageNum, (Page *&)node);
		while (node->level == 0) {
			const int idx = _search(node, lowValInt);
			const PageId oldPageNum = pageNum;
			pageNum = node->pageNoArray[idx];
			bufMgr->unPinPage(file, oldPageNum, false);
			bufMgr->readPage(file, pageNum, (Page *&)node);
		}
		{
			const int idx = _search(node, lowValInt);
			const PageId oldPageNum = pageNum;
			pageNum = node->pageNoArray[idx];
			bufMgr->unPinPage(file, oldPageNum, false);
			bufMgr->readPage(file, pageNum, (Page *&)leafNode);
		}

		bool found = false;

		for (int idx = 0; idx < INTARRAYLEAFSIZE && leafNode->keyArray[idx] != EMPTY_SLOT; idx++) {
			const int key = leafNode->keyArray[idx];
			if (_satisfies(lowValInt, lowOpParm, highValInt, highOpParm, key)) {
				found = true;
				currentPageNum = pageNum;	
				currentPageData = (Page *)leafNode;
				nextEntry = idx;
				scanExecuting = true;
				break;
			}
		}

		if (!found) {
			scanExecuting = false;
			bufMgr->unPinPage(file, pageNum, false);
			throw NoSuchKeyFoundException();
		}
	}

	// -----------------------------------------------------------------------------
	// BTreeIndex::scanNext
	// -----------------------------------------------------------------------------

	const void BTreeIndex::scanNext(RecordId& outRid) 
	{
		if (!scanExecuting) {
			throw ScanNotInitializedException();	
		}

		LeafNodeInt *node = (LeafNodeInt *)currentPageData;


		{
/*
			std::cout << "leaf node" << std::endl;
			std::cout << "next Entry: " << nextEntry << std::endl;
			std::cout << "next sib:" << node->rightSibPageNo << std::endl;
			for (int i = 0; i < INTARRAYLEAFSIZE; i++) {
				std::cout << node->keyArray[i] << ", ";
			}
			std::cout << std::endl;
*/
		}

		if (nextEntry == INTARRAYLEAFSIZE || node->keyArray[nextEntry] == EMPTY_SLOT) {
			if (node->rightSibPageNo == (PageId)EMPTY_SLOT) {
				throw IndexScanCompletedException();
			}
			const PageId oldPageNum = currentPageNum;
			currentPageNum = node->rightSibPageNo;
			nextEntry = 0;
			bufMgr->unPinPage(file, oldPageNum, false);
			bufMgr->readPage(file, currentPageNum, (Page *&)currentPageData);
			node = (LeafNodeInt *)currentPageData;
		}
		
		const int key = node->keyArray[nextEntry];
		
//		std::cout << key << ", "<< nextEntry << std::endl;
		if (_satisfies(lowValInt, lowOp, highValInt, highOp, key)) {
			outRid = node->ridArray[nextEntry];
			nextEntry++;
		} else {
			throw IndexScanCompletedException();
		}
	}

	// -----------------------------------------------------------------------------
	// BTreeIndex::endScan
	// -----------------------------------------------------------------------------
	//
	const void BTreeIndex::endScan() 
	{
		if ( !scanExecuting) {
			throw ScanNotInitializedException();	
		}
		bufMgr->unPinPage(file, currentPageNum, false);

		scanExecuting = false;
		currentPageData = nullptr;
		currentPageNum = -1;
		nextEntry = -1;
	}

	const bool BTreeIndex::_satisfies(int lowVal, const Operator lowOp, int highVal, const Operator highOp, int val)
	{
		if (lowOp == GTE && highOp == LTE) {
			return val >= lowVal && val <= highVal;
		} else if (lowOp == GT && highOp == LTE) {
			return val > lowVal && val <= highVal;
		} else if (lowOp == GTE && highOp == LT) {
			return val >= lowVal && val < highVal;
		} else {
			return val > lowVal && val < highVal;
		}
	}

	const void BTreeIndex::_checkScanParams(int lowVal, const Operator lowOp, int highVal, const Operator highOp) 
	{
		if (lowOp != GT && lowOp != GTE) {
			throw BadOpcodesException();
		}
		if (highOp != LT && highOp != LTE) {
			throw BadOpcodesException();
		}
		if (lowVal > highVal) {
			throw BadScanrangeException();
		}
	}

	const int BTreeIndex::_search(NonLeafNodeInt *node, int val) 
	{
		int idx = 0;
		for (;idx < INTARRAYNONLEAFSIZE && node->keyArray[idx] < val && node->keyArray[idx] != EMPTY_SLOT; idx++);
		return idx;
	}

	const void BTreeIndex::_assertLeafInternalConsistency(LeafNodeInt *node)
	{
		for (int i = 0; i < INTARRAYLEAFSIZE; i++)
		{
			if (node->keyArray[i] == EMPTY_SLOT) {
				assert(node->ridArray[i].page_number == (PageId)EMPTY_SLOT);
				assert(node->ridArray[i].slot_number == (SlotId)EMPTY_SLOT);
			} else {
				assert(node->ridArray[i].page_number != (PageId)EMPTY_SLOT);
				assert(node->ridArray[i].slot_number != (SlotId)EMPTY_SLOT);
			}
		}
	}

	const void BTreeIndex::_assertNonLeafInternalConsistency(NonLeafNodeInt *node)
	{
		assert(node->level == 0 || node->level == 1);
		for (int i = 0; i < INTARRAYNONLEAFSIZE; i++)
		{
			if (node->keyArray[i] == EMPTY_SLOT) {
				assert(node->pageNoArray[i+1] == (PageId)EMPTY_SLOT);
			} else {
				assert(node->pageNoArray[i+1] != (PageId)EMPTY_SLOT);
			}
		}
	}

	PageIDPair BTreeIndex::_newLeafNode()
	{
		PageId pid;
		Page *page;
		bufMgr->allocPage(file, pid, page);

		LeafNodeInt *node = reinterpret_cast<LeafNodeInt *>(page);

		node -> rightSibPageNo = EMPTY_SLOT;
		for (int i = 0; i < INTARRAYLEAFSIZE; i++) {
			node->keyArray[i] = EMPTY_SLOT;
			node->ridArray[i].page_number = (PageId)EMPTY_SLOT;
			node->ridArray[i].slot_number = (SlotId)EMPTY_SLOT;
		}
		
		PageIDPair pair;
		pair.set(page, pid);
		return pair;
	}

	PageIDPair BTreeIndex::_newNonLeafNode()
	{
		PageId pid;
		Page *page;
		bufMgr->allocPage(file, pid, page);

		NonLeafNodeInt *node = reinterpret_cast<NonLeafNodeInt *>(page);

		node -> level = 0;
		for (int i = 0; i < INTARRAYNONLEAFSIZE; i++) {
			node->keyArray[i] = (PageId)EMPTY_SLOT;
			node->pageNoArray[i] = (PageId)EMPTY_SLOT;
		}
		node->pageNoArray[INTARRAYNONLEAFSIZE] = (PageId)EMPTY_SLOT;

		PageIDPair pair;
		pair.set(page, pid);
		return pair;
	}

	const void BTreeIndex::_internalTest()
	{
		_testValidateMetaPage();
		_testGetRIDKeyPair();
		_testNewPageIsConsistent();
		_testLeafInsertEntry();
		_testNonLeafInsertEntry();
		_testLeafSplitInsertEntry();
		_testNonLeafSplitInsertEntry();
	}

	const void BTreeIndex::_testValidateMetaPage()
	{
		IndexMetaInfo meta;
		meta.relationName[0] = 'p';
		meta.relationName[1] = 0;	
		meta.attrByteOffset = 10;
		meta.attrType = Datatype::INTEGER;
		meta.rootPageNo = 4;	

		std::string name = "p";
		assert(_validateMetaPage(&meta, name, 10, Datatype::INTEGER));
		assert(!_validateMetaPage(&meta, name, 9, Datatype::INTEGER));
		std::cout << "test validate meta page passes" << std::endl;
	}

	const void BTreeIndex::_testGetRIDKeyPair()
	{
		{
			char key[4];
			*((int *)key) = 10;
			std::string rec(reinterpret_cast<char*>(key), 4);

			RecordId rid;
			rid.page_number = 2;
			rid.slot_number = 4;
			RIDKeyPair<int> res = _getRIDKeyPairFromRecord(rec, rid, 0);
			assert(res.rid == rid);
			assert(res.key == 10);
		}

		{
			char key2[4];
			*((int *)key2) = 10;
			std::string rec2 = "1234";
			std::string tmp(reinterpret_cast<char*>(key2), 4);
			rec2 = rec2 + tmp;
			assert(rec2.length() == 8);

			RecordId rid2;
			rid2.page_number = 2;
			rid2.slot_number = 4;
			RIDKeyPair<int> res2 = _getRIDKeyPairFromRecord(rec2, rid2, 4);
			assert(res2.rid == rid2);
			assert(res2.key == 10);
		}

		std::cout << "test get rid key pair passes" << std::endl;
	}

	const void BTreeIndex::_testNewPageIsConsistent()
	{
		{
			PageIDPair pair = _newLeafNode();
			_assertLeafInternalConsistency(reinterpret_cast<LeafNodeInt *>(pair.page));
			bufMgr->unPinPage(file, pair.pageNo, false);
		}
		{
			PageIDPair pair = _newNonLeafNode();
			_assertNonLeafInternalConsistency(reinterpret_cast<NonLeafNodeInt *>(pair.page));
			bufMgr->unPinPage(file, pair.pageNo, false);
		}
		std::cout << "test new node is consistent passed" << std::endl;
	}

	const void BTreeIndex::_testLeafInsertEntry()
	{
		{
			PageIDPair pair = _newLeafNode();
			LeafNodeInt *node = reinterpret_cast<LeafNodeInt *>(pair.page);
			int k1 = 1, k2 = 2, k3 = 3, k4 = 4;
			RecordId r1 = {1,1}, r2 = {1,2}, r3 = {1,3}, r4 = {1,4};
			_leafInsertEntry(node, &k2, r2);
			_leafInsertEntry(node, &k1, r1);
			_leafInsertEntry(node, &k4, r4);
			_leafInsertEntry(node, &k3, r3);
			for (int i = 0; i < 4; i++) {
				assert(node->keyArray[i] == i+1);
				assert(node->ridArray[i].slot_number == i+1);
			}
			bufMgr->unPinPage(file, pair.pageNo, false);
		}
		std::cout << "test lead insert entry passes" << std::endl;
	}

	const void BTreeIndex::_testNonLeafInsertEntry()
	{
		{
			PageIDPair pair = _newNonLeafNode();
			NonLeafNodeInt *node = reinterpret_cast<NonLeafNodeInt *>(pair.page);
			node->pageNoArray[0] = 987;
			PageKeyPair<int> p1, p2, p3, p4;
			p1.set(1, 1);
			p2.set(2, 2);
			p3.set(3, 3);
			p4.set(4, 4);
			_nonLeafInsertEntry(node, p2);
			_nonLeafInsertEntry(node, p1);
			_nonLeafInsertEntry(node, p4);
			_nonLeafInsertEntry(node, p3);
			assert(node->pageNoArray[0] == 987);
			for (int i = 0; i < 4; i++) {
				assert(node->keyArray[i] == i+1);
				assert(node->pageNoArray[i+1] == PageId(i+1));
			}
			bufMgr->unPinPage(file, pair.pageNo, false);
		}
		std::cout << "test non leaf insert entry passes" << std::endl;
	}

	const void BTreeIndex::_testLeafSplitInsertEntry()	{
		{
			PageIDPair pair = _newLeafNode();
			LeafNodeInt *node = reinterpret_cast<LeafNodeInt *>(pair.page);
			node->rightSibPageNo = 756;
			for (int i = 0; i < INTARRAYLEAFSIZE; i++) {
				node->keyArray[i] = i + 1;
				node->ridArray[i] = {1, (SlotId) (i + 1)};
			}
			int k1 = 1;
			RecordId r1 = {1,1};
			PageKeyPair<int> newNode = _leafSplitInsertEntry(node, &k1, r1);
			assert(node->keyArray[0] == 1);
			for (int i = 1; i < INTARRAYLEAFSIZE / 2;  i++) {
				assert(node->keyArray[i] == i);
				assert(node->ridArray[i].slot_number == i);
			}

			for (int i = INTARRAYLEAFSIZE / 2;  i < INTARRAYLEAFSIZE; i++) {
				assert(node->keyArray[i] == EMPTY_SLOT);
				assert(node->ridArray[i].slot_number == (SlotId) EMPTY_SLOT);
			}
			Page* page;
			bufMgr->readPage(file, newNode.pageNo, page);
			LeafNodeInt *newNode1 = reinterpret_cast<LeafNodeInt *>(page);
			for (int i = 0; i <= INTARRAYLEAFSIZE / 2;  i++) {
				assert(newNode1->keyArray[i] == INTARRAYLEAFSIZE / 2 + i);
				assert(newNode1->ridArray[i].slot_number == INTARRAYLEAFSIZE / 2 + i);
			}

			for (int i = INTARRAYLEAFSIZE / 2 + 1;  i < INTARRAYLEAFSIZE; i++) {
				assert(newNode1->keyArray[i] == EMPTY_SLOT);
				assert(newNode1->ridArray[i].slot_number == (SlotId) EMPTY_SLOT);
			}
		
			assert(node->rightSibPageNo == newNode.pageNo);
			assert(newNode1->rightSibPageNo == 756);
			bufMgr->unPinPage(file, pair.pageNo, false);
			bufMgr->unPinPage(file, newNode.pageNo, false);
		}

		{
			PageIDPair pair = _newLeafNode();
			LeafNodeInt *node = reinterpret_cast<LeafNodeInt *>(pair.page);
			node->rightSibPageNo = 756;
			for (int i = 0; i < INTARRAYLEAFSIZE; i++) {
				node->keyArray[i] = i + 1;
				node->ridArray[i] = {1, (SlotId) (i + 1)};
			}
			int k1 = INTARRAYLEAFSIZE / 2 - 1;;
			RecordId r1 = {(PageId) k1, (SlotId) k1};
			PageKeyPair<int> newNode = _leafSplitInsertEntry(node, &k1, r1);
			assert(node->keyArray[INTARRAYLEAFSIZE / 2 - 1] == k1);
			for (int i = 0; i < INTARRAYLEAFSIZE / 2 - 2;  i++) {
				assert(node->keyArray[i] == i + 1);
				assert(node->ridArray[i].slot_number == i + 1);
			}

			for (int i = INTARRAYLEAFSIZE / 2;  i < INTARRAYLEAFSIZE; i++) {
				assert(node->keyArray[i] == EMPTY_SLOT);
				assert(node->ridArray[i].slot_number == (SlotId) EMPTY_SLOT);
			}
			Page* page;
			bufMgr->readPage(file, newNode.pageNo, page);
			LeafNodeInt *newNode1 = reinterpret_cast<LeafNodeInt *>(page);
			for (int i = 0; i <= INTARRAYLEAFSIZE / 2;  i++) {
				assert(newNode1->keyArray[i] == INTARRAYLEAFSIZE / 2 + i);
				assert(newNode1->ridArray[i].slot_number == INTARRAYLEAFSIZE / 2 + i);
			}

			for (int i = INTARRAYLEAFSIZE / 2 + 1;  i < INTARRAYLEAFSIZE; i++) {
				assert(newNode1->keyArray[i] == EMPTY_SLOT);
				assert(newNode1->ridArray[i].slot_number == (SlotId) EMPTY_SLOT);
			}
		
			assert(node->rightSibPageNo == newNode.pageNo);
			assert(newNode1->rightSibPageNo == 756);
			bufMgr->unPinPage(file, pair.pageNo, false);
			bufMgr->unPinPage(file, newNode.pageNo, false);
		}

		{
			PageIDPair pair = _newLeafNode();
			LeafNodeInt *node = reinterpret_cast<LeafNodeInt *>(pair.page);
			node->rightSibPageNo = 756;
			for (int i = 0; i < INTARRAYLEAFSIZE; i++) {
				node->keyArray[i] = i + 1;
				node->ridArray[i] = {1, (SlotId) (i + 1)};
			}
			int k1 = INTARRAYLEAFSIZE / 2 + 1;
			RecordId r1 = {(PageId) k1, (SlotId) k1};
			PageKeyPair<int> newNode = _leafSplitInsertEntry(node, &k1, r1);
			for (int i = 0; i < INTARRAYLEAFSIZE / 2;  i++) {
				assert(node->keyArray[i] == i + 1);
				assert(node->ridArray[i].slot_number == i + 1);
			}

			for (int i = INTARRAYLEAFSIZE / 2;  i < INTARRAYLEAFSIZE; i++) {
				assert(node->keyArray[i] == EMPTY_SLOT);
				assert(node->ridArray[i].slot_number == (SlotId) EMPTY_SLOT);
			}
			Page* page;
			bufMgr->readPage(file, newNode.pageNo, page);
			LeafNodeInt *newNode1 = reinterpret_cast<LeafNodeInt *>(page);	
			assert(newNode1->keyArray[0] == k1);
			for (int i = 1; i <= INTARRAYLEAFSIZE / 2;  i++) {
				assert(newNode1->keyArray[i] == INTARRAYLEAFSIZE / 2 + i);
				assert(newNode1->ridArray[i].slot_number == INTARRAYLEAFSIZE / 2 + i);
			}

			for (int i = INTARRAYLEAFSIZE / 2 + 1;  i < INTARRAYLEAFSIZE; i++) {
				assert(newNode1->keyArray[i] == EMPTY_SLOT);
				assert(newNode1->ridArray[i].slot_number == (SlotId) EMPTY_SLOT);
			}
		
			assert(node->rightSibPageNo == newNode.pageNo);
			assert(newNode1->rightSibPageNo == 756);
			bufMgr->unPinPage(file, pair.pageNo, false);
			bufMgr->unPinPage(file, newNode.pageNo, false);
		}
		{
			PageIDPair pair = _newLeafNode();
			LeafNodeInt *node = reinterpret_cast<LeafNodeInt *>(pair.page);
			node->rightSibPageNo = 756;
			for (int i = 0; i < INTARRAYLEAFSIZE; i++) {
				node->keyArray[i] = i + 1;
				node->ridArray[i] = {1, (SlotId) (i + 1)};
			}
			int k1 = 10000;
			RecordId r1 = {(PageId) k1, (SlotId) k1};
			PageKeyPair<int> newNode = _leafSplitInsertEntry(node, &k1, r1);
			for (int i = 0; i < INTARRAYLEAFSIZE / 2;  i++) {
				assert(node->keyArray[i] == i + 1);
				assert(node->ridArray[i].slot_number == i + 1);
			}

			for (int i = INTARRAYLEAFSIZE / 2;  i < INTARRAYLEAFSIZE; i++) {
				assert(node->keyArray[i] == EMPTY_SLOT);
				assert(node->ridArray[i].slot_number == (SlotId) EMPTY_SLOT);
			}
			Page* page;
			bufMgr->readPage(file, newNode.pageNo, page);
			LeafNodeInt *newNode1 = reinterpret_cast<LeafNodeInt *>(page);	
			assert(newNode1->keyArray[INTARRAYLEAFSIZE / 2] == k1);
			for (int i = 0; i < INTARRAYLEAFSIZE / 2;  i++) {
				assert(newNode1->keyArray[i] == INTARRAYLEAFSIZE / 2 + i + 1);
				assert(newNode1->ridArray[i].slot_number == INTARRAYLEAFSIZE / 2 + i + 1);
			}

			for (int i = INTARRAYLEAFSIZE / 2 + 1;  i < INTARRAYLEAFSIZE; i++) {
				assert(newNode1->keyArray[i] == EMPTY_SLOT);
				assert(newNode1->ridArray[i].slot_number == (SlotId) EMPTY_SLOT);
			}
		
			assert(node->rightSibPageNo == newNode.pageNo);
			assert(newNode1->rightSibPageNo == 756);
			bufMgr->unPinPage(file, pair.pageNo, false);
			bufMgr->unPinPage(file, newNode.pageNo, false);
		}
		std::cout << "test leaf split insert entry passes" << std::endl;
	}
	

	const void BTreeIndex::_testNonLeafSplitInsertEntry()
	{
		{
			PageIDPair pair = _newNonLeafNode();
			NonLeafNodeInt *node = reinterpret_cast<NonLeafNodeInt *>(pair.page);
			for (int i = 0; i < INTARRAYNONLEAFSIZE; i++) {
				node->keyArray[i] = i + 1;
				node->pageNoArray[i + 1] = i + 1;
			}
			node->pageNoArray[0] = 0;

			int k1 = 1;
			const int half = (INTARRAYNONLEAFSIZE + 1) / 2;
			PageKeyPair<int> pairIn;
			pairIn.set(k1, k1);

			PageKeyPair<int> newNodePair = _nonLeafSplitInsertEntry(node, pairIn);
			assert(node->pageNoArray[0] == (PageId)(0));

			assert(node->keyArray[0] == 1);
			assert(node->pageNoArray[1] == (PageId)(1));
			
			for (int i = 1; i < half;  i++) {
				assert(node->keyArray[i] == i);
				//std::cout << "expect: " << i << "actual: " << node->pageNoArray[i] << std::endl;
				assert(node->pageNoArray[i+1] == (PageId)(i));
			}
			assert(node->pageNoArray[half] == (PageId)(half - 1));
			for (int i = half;  i < INTARRAYNONLEAFSIZE; i++) {
				assert(node->keyArray[i] == EMPTY_SLOT);
				assert(node->pageNoArray[i+1] == (PageId)(EMPTY_SLOT));
			}

			Page* page;
			bufMgr->readPage(file, newNodePair.pageNo, page);
			NonLeafNodeInt *newNode = reinterpret_cast<NonLeafNodeInt *>(page);	
			assert(newNode->pageNoArray[0] == PageId(half));
			for (int i = 0; i < half - 1;  i++) {
				assert(newNode->keyArray[i] == half + i + 1);
				assert(newNode->pageNoArray[i + 1] == PageId(half + i + 1));
			}
			assert(newNode->pageNoArray[half - 1] == (PageId)(INTARRAYNONLEAFSIZE));
			for (int i = half;  i < INTARRAYNONLEAFSIZE; i++) {
				assert(newNode->keyArray[i] == EMPTY_SLOT);
				assert(newNode->pageNoArray[i + 1] == (PageId)EMPTY_SLOT);
			}

			assert(newNodePair.key == half);

			bufMgr->unPinPage(file, pair.pageNo, false);
			bufMgr->unPinPage(file, newNodePair.pageNo, false);
		}
		{
			const int half = (INTARRAYNONLEAFSIZE + 1) / 2;

			PageIDPair pair = _newNonLeafNode();
			NonLeafNodeInt *node = reinterpret_cast<NonLeafNodeInt *>(pair.page);
			for (int i = 0; i < half; i++) {
				node->keyArray[i] = i + 1;
				node->pageNoArray[i + 1] = i + 1;
			}
			for (int i = half; i < INTARRAYNONLEAFSIZE; i++) {
				node->keyArray[i] = i + 2;
				node->pageNoArray[i + 1] = i + 2;
			}
			node->pageNoArray[0] = 0;

			int k1 = half + 1;
			PageKeyPair<int> pairIn;
			pairIn.set(k1, k1);

			PageKeyPair<int> newNodePair = _nonLeafSplitInsertEntry(node, pairIn);
			assert(node->pageNoArray[0] == (PageId)(0));
			for (int i = 0; i < half;  i++) {
				assert(node->keyArray[i] == i + 1);
				//std::cout << "expect: " << i << "actual: " << node->pageNoArray[i] << std::endl;
				assert(node->pageNoArray[i+1] == (PageId)(i + 1));
			}
			assert(node->pageNoArray[half] == (PageId)(half));
			for (int i = half;  i < INTARRAYNONLEAFSIZE; i++) {
				assert(node->keyArray[i] == EMPTY_SLOT);
				assert(node->pageNoArray[i+1] == (PageId)(EMPTY_SLOT));
			}

			Page* page;
			bufMgr->readPage(file, newNodePair.pageNo, page);
			NonLeafNodeInt *newNode = reinterpret_cast<NonLeafNodeInt *>(page);	
			assert(newNode->pageNoArray[0] == PageId(k1));
			for (int i = 0; i < half - 1;  i++) {
				assert(newNode->keyArray[i] == half + i + 2);
				assert(newNode->pageNoArray[i + 1] == PageId(half + i + 2));
			}
			for (int i = half - 1;  i < INTARRAYNONLEAFSIZE; i++) {
				assert(newNode->keyArray[i] == EMPTY_SLOT);
				assert(newNode->pageNoArray[i + 1] == (PageId)EMPTY_SLOT);
			}

			assert(newNodePair.key == k1);

			bufMgr->unPinPage(file, pair.pageNo, false);
			bufMgr->unPinPage(file, newNodePair.pageNo, false);
		}		
		{
			PageIDPair pair = _newNonLeafNode();
			NonLeafNodeInt *node = reinterpret_cast<NonLeafNodeInt *>(pair.page);
			for (int i = 0; i < INTARRAYNONLEAFSIZE; i++) {
				node->keyArray[i] = i + 1;
				node->pageNoArray[i + 1] = i + 1;
			}
			node->pageNoArray[0] = 0;

			int k1 = 9999;
			const int half = (INTARRAYNONLEAFSIZE + 1) / 2;
			PageKeyPair<int> pairIn;
			pairIn.set(k1, k1);

			PageKeyPair<int> newNodePair = _nonLeafSplitInsertEntry(node, pairIn);
			assert(node->pageNoArray[0] == (PageId)(0));

			for (int i = 0; i < half;  i++) {
				assert(node->keyArray[i] == i + 1);
				assert(node->pageNoArray[i+1] == (PageId)(i + 1));
			}
			for (int i = half;  i < INTARRAYNONLEAFSIZE; i++) {
				assert(node->keyArray[i] == EMPTY_SLOT);
				assert(node->pageNoArray[i+1] == (PageId)(EMPTY_SLOT));
			}

			Page* page;
			bufMgr->readPage(file, newNodePair.pageNo, page);
			NonLeafNodeInt *newNode = reinterpret_cast<NonLeafNodeInt *>(page);	
			assert(newNode->pageNoArray[0] == PageId(half + 1));
			for (int i = 0; i < half - 2;  i++) {
				assert(newNode->keyArray[i] == half + i + 2);
				assert(newNode->pageNoArray[i + 1] == PageId(half + i + 2));
			}
			assert(newNode->keyArray[half - 2] == 9999);
			assert(newNode->pageNoArray[half - 1] == PageId(9999));
			for (int i = half - 1;  i < INTARRAYNONLEAFSIZE; i++) {
				assert(newNode->keyArray[i] == EMPTY_SLOT);
				assert(newNode->pageNoArray[i + 1] == PageId(EMPTY_SLOT));
			}

			assert(newNodePair.key == half + 1);

			bufMgr->unPinPage(file, pair.pageNo, false);
			bufMgr->unPinPage(file, newNodePair.pageNo, false);
		}
		std::cout << "test non leaf split insert entry passes" << std::endl;
	}

}
