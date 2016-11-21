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

#define EMPTY_SLOT (-1)

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
		}

		// If the index file does not exist, then a new file is created
		catch (FileNotFoundException e)
		{
			// If the file does not exist, a new index file is created
			file = new BlobFile(outIndexName, true);
			std::cout << "A new index file is created";// Test output message    

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
					// std::cout << "My record: " << record << std::endl;
					// std::cout << outRid.page_number << "   " << outRid.slot_number << "\n";

				}
			}
			catch (EndOfFileException e)
			{
				std::cout << "Reach the end of the file." << std::endl;
			}
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
	}

	// -----------------------------------------------------------------------------
	// BTreeIndex::insertEntry
	// -----------------------------------------------------------------------------

	const void BTreeIndex::insertEntry(const void *key, const RecordId rid) 
	{

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

	// -----------------------------------------------------------------------------
	// BTreeIndex::startScan
	// -----------------------------------------------------------------------------

	const void BTreeIndex::startScan(const void* lowValParm,
			const Operator lowOpParm,
			const void* highValParm,
			const Operator highOpParm)
	{

	}

	// -----------------------------------------------------------------------------
	// BTreeIndex::scanNext
	// -----------------------------------------------------------------------------

	const void BTreeIndex::scanNext(RecordId& outRid) 
	{

	}

	// -----------------------------------------------------------------------------
	// BTreeIndex::endScan
	// -----------------------------------------------------------------------------
	//
	const void BTreeIndex::endScan() 
	{

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

	LeafNodeInt *BTreeIndex::_newLeafNode()
	{
		LeafNodeInt *node = (LeafNodeInt *)malloc(sizeof(LeafNodeInt));
		node -> rightSibPageNo = EMPTY_SLOT;
		for (int i = 0; i < INTARRAYLEAFSIZE; i++) {
			node->keyArray[i] = EMPTY_SLOT;
			node->ridArray[i].page_number = (PageId)EMPTY_SLOT;
			node->ridArray[i].slot_number = (SlotId)EMPTY_SLOT;
		}
		return node;
	}

	NonLeafNodeInt *BTreeIndex::_newNonLeafNode()
	{
		NonLeafNodeInt *node = (NonLeafNodeInt *)malloc(sizeof(NonLeafNodeInt));
		node -> level = 0;
		for (int i = 0; i < INTARRAYNONLEAFSIZE; i++) {
			node->keyArray[i] = (PageId)EMPTY_SLOT;
			node->pageNoArray[i] = (PageId)EMPTY_SLOT;
		}
		node->pageNoArray[INTARRAYNONLEAFSIZE] = (PageId)EMPTY_SLOT;
		return node;
	}

	const void BTreeIndex::_internalTest()
	{
		_testValidateMetaPage();
		_testGetRIDKeyPair();
		_testNewPageIsConsistent();
		_testLeafInsertEntry();
		_testNonLeafInsertEntry();
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
			LeafNodeInt *ln = _newLeafNode();
			_assertLeafInternalConsistency(ln);
			free(ln);
		}
		{
			NonLeafNodeInt *n = _newNonLeafNode();
			_assertNonLeafInternalConsistency(n);
			free(n);
		}
		std::cout << "test new node is consistent passed" << std::endl;
	}

	const void BTreeIndex::_testLeafInsertEntry()
	{
		{
			LeafNodeInt *node = _newLeafNode();
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
			free(node);
		}
		std::cout << "test lead insert entry passes" << std::endl;
	}

	const void BTreeIndex::_testNonLeafInsertEntry()
	{
		{
			NonLeafNodeInt *node = _newNonLeafNode();
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
			free(node);
		}
		std::cout << "test non leaf insert entry passes" << std::endl;
	}
}
