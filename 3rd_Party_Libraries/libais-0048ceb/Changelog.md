Changelog and release history for libais
========================================

0.15 - 2015-06-16
-----------------

* Added libais namespace
* Started using clang-format
* Fix Ais18 carrier sense (CS) messages
* Use AisBitset for bit decoding
* Fix AIS_STATUS_STRINGS
* Use std::array for fixed sized arrays
* Update DAC enum with new values
* Added ais/compatibility/gpsd.py
* Added ais/stream (deprecated - please do not use)
* Closed github issues 27, 31, 34, 36, 38, 42, 45, 46, 47, 51, 52, 53, 56, and 60
* Added ais/nmea.py, ais/tag_block.py, and ais/util_test.py
* Rearranged the tree structure to have src/... for c++ and {ais,test}/ for python
* Use more c++11: NULL -&gt; nullptr
* Add 8:367:22 Area Notice
* python setup.py test mostly works
* At least one test for all top levels except msg 6.
* Add Travis CI testing

0.14 - 2014-04-22
-----------------

* Switch license from LGPL 3+ to Apache 2.0
* Msg 8:366:22 bit count check fixed
* Fixed spare bit calculation in msg 8
* Added DAC country codes and FI message ids enums
* Fixed error tracking in msg 6 and 8 sub-messages

0.13 - 2012-11-18
-----------------

* Switch to the [Google C++ style guide](http://code.google.com/p/google-styleguide/) and [`cpplint.py`](http://google-styleguide.googlecode.com/svn/trunk/cpplint/)
* Lots of small bugs found in the code review process
* Message constructors now start as status = AIS_UNINITIALIZED and set to AIS_OK when decoding is done
* Switched to using initializers to call parent initializers in C++
* Removed a lot of duplicate and unused code. AIS_ERR_WRONG_MSG_TYPE removed.  More still needs to be removed.
* Rewrote the C++ NMEA parsing functions: GetBody, GetNthField, GetPad and Split.
* Switch to using asserts for checks that imply coding errors within the library.


0.12 - 2012-11-05
-----------------

* Fix bit count bugs in 8_1_14, 8_1_15, 8_1_27
* Rewrote nth_field.  Added split and delimiters
* Folded in 366 header to ais.h
* Lots and lots of style cleanup
* Use std::foo, remove std:: from code
* Documentation for Msg 17
* Testing of Msg 20


0.11 - 2012-10-29
-----------------

* New release because a large binary went out in 0.10


0.10 - 2012-10-29
-----------------

* Add a test directory and test of all top level msgs except 20 in python
* Begin cleanup of test_libais.cpp
* Almost all decoders require pad bits now
* linted - lots of formatting changes
* Start implmenting Msg 17 GNSS differential corrections
* Convert FIX to TODO and put (schwehr) after each to assign them to myself.
* remove bool casting of bitset[offset]
* Implemented 8 1 26
* Clean up c++ logical oprators.  Do not use and, or, and not
* Message 24 needed pad.  Fixed
* Removed print()
* remove almost all cout/cerr that were not in print()
* Remove lots of dead code
* Pass pad into ais.decode in python, but handle without


0.9 - 2012-10-19
-----------------

* ais.decode now requires the pad bits in python
* Added RIS 8_200_{10,23,24,40,55}
* Implmented the rest of Circ 236 BBM
* Implmented all Circ 289 messages except ABM route and BBM env sensors
* Implemented AIS messages 6, 9, 10, 12, 15-17, 20-22, 23, 25-27. Still payload work to do.
* Imported rolker's CMakeList.txt


0.8 - 2012-05-12
-----------------

* MANIFEST.in now has VERSION
* Implemented AIS messages 7, 13, 14, 18, 19, 24


0.7 - 2012-04-30
-----------------

* Added MANIFEST.in
* setup.py compliant with pypi
* Added AIS area/zone messages 8:1:22 and 8:366:22


0.6 - 2010-06-17
-----------------

* ais21.cpp: new file - AtoN status
* ais.h: fix CHECKPOINT for emacs 23
* ais.h: Proper inheritance from AisMsg of message_id, repead_id and mmsi
* ais.h: started trying to define 8_366_34 - zone msg
* nagios_pg_ais.py: new file - monitor db statios with nagios ssh or snmp


0.5 - 2010-05-25
-----------------

* Still a lot of untested/unimplemented messages
* docs: Include ESR's AIVDM.txt with permission
* docs: MID / DAC / MMSI prefixes now listed in mid.csv
* docs: dac/fi list
* docs: More notes for message designers
* Added msg 8 - 1:11 - IMO Met/Hydro
* Added AIS msg 9 - SAR Position
* nais2py.py: Added try except wrapper on ProcessingThread.  Also try to track one off error found on call with timestamp converting to float
* send_data.py: new file for testing
* nais2py.py: LineQueue now has a custom drop handler for too many lines waiting.  Could be better.  
* nais2py.py: Added threaded network interface.  Seeing the network side overwhelm the processing thread
* nais2py.py: Added response_class handling to VesselNames.  Can be preloaded.  Allows periodic name dump
* nais2py.py: Added ENABLE_DB flag to try runs without database execute commands.  Faster debugging
* vessels.csv: new file - example preloading of vessel names as response ships


0.4 - 2010-05-11
-----------------

* nais2py.py: Started providing a command line interface
* nais2py.py: Added PositionCache class
* nais2py.py: Added distance calculation code.  
* nais2py.py: Changed the database table names and structure.  Now vessel_name and vessel_pos
* ais_lut.py: new file with lookup tables to make ais wire codes human readable.


0.3 - 2010-05-10
-----------------

* ais.c: added check_error_messages to make sure they are not out of sync
* -D_GLIBCXX_DEBUG appears broken in GCC 4.[0-2] so do not use
* ais_decode_normed.cpp: temporary C++ side decoding of normed AIVDM messages
* nais2pg: added vesselname class to manage updates to postgresql
* Added message 24
* Fixed python reference counting.  Added XXSafeSetItem functions


0.2 - 2010-05-06
-----------------

* Added C++ error handling to classes via AIS_STATUS
* C++ message now inherit from AisMsg and need to call init() in constructor
* Added C++ messages 7_13, 14, 18, and 19
* aivdm_to_bits now has error checking
* ais123.cpp renames to ais1_2_3.cpp
* Switched to unicode in ais_py.cpp to support Python 3
* ais_py.cpp has strange INIT to handle Python 2 and 3
* nais2pg.py is starting to implement a new AIS feed to database daemon
* Redid my old USCG regex to have better names with lower_lower style
* LineQueue should now support reading through a socket, which I got wrong before
* Total redo of the Normalization queue to be much lower overhead and to call the regex only once per line received.  Only keep they body of all but the last message in a sequence.
* test_libais.cpp is not really a test framework, but it does try out the pure C++ world.


0.1 - 2010-05-03
-----------------

* Able to decode messages 1-5 from python
* Still a lot of work left to do!
