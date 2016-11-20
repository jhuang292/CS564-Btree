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

	const void BTreeIndex::_internalTest()
	{
		_testValidateMetaPage();
		_testGetRIDKeyPair();
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

}
