============
Introduction
============

.. image:: https://travis-ci.org/schwehr/libais.svg?branch=master
    :target: https://travis-ci.org/schwehr/libais

.. image:: https://scan.coverity.com/projects/5519/badge.svg
    :target: https://scan.coverity.com/projects/5519

.. image:: https://codeclimate.com/github/schwehr/libais/badges/gpa.svg
    :target: https://codeclimate.com/github/schwehr/libais

.. image:: https://badge.fury.io/py/libais.svg
    :target: http://badge.fury.io/py/libais

Library for decoding maritime Automatic Identification System messages.

See Also
========

`Automatic Identification System <http://en.wikipedia.org/wiki/Automatic_Identification_System>`_

Other open source AIS projects:

- `GPSd <http://en.wikipedia.org/wiki/Gpsd>`_
- `FreeAIS.org <http://www.freeais.org/>`_ - does this actually decode AIS?
- `AisLib <http://github.com/DaMSA/AisLib>`_
- `noaadata <http://github.com/schwehr/noaadata>`_
- `ais-areanotice <https://github.com/schwehr/ais-areanotice-py>`_
- `OpenCPN <https://github.com/OpenCPN/OpenCPN>`_
- `aisparser <https://github.com/bcl/aisparser>`_

Building
========

Building with legacy Makefile
-----------------------------

.. code-block:: console

    $ make -f Makefile-custom test

Building with Python
--------------------

.. code-block:: console

    $ python setup.py build
    $ python setup.py install

Building with CMake
-------------------

.. code-block:: console

    $ cmake .
    $ make

Usage
=====

There are two interfaces to libais, one high-level iterator based one
and a low-level fast C++ only one. The iterator based interface is
accessed the following way:

.. code-block:: python

    import ais.stream
    with open("file.nmea") as f:
        for msg in ais.stream.decode(f):
            print msg

To use the low-level C++ interface directly, you need to handle multi-line messages and padding yourself:

.. code-block:: python

    import ais
    ais.decode('15PIIv7P00D5i9HNn2Q3G?wB0t0I', 0)
    ais.decode('402u=TiuaA000r5UJ`H4`?7000S:', 0)
    ais.decode('55NBjP01mtGIL@CW;SM<D60P5Ld000000000000P0`<3557l0<50@kk@K5h@00000000000', 0)

There is also support for converting parsed messages to the structure
output by GPSD / gpsdecode. For full compatibility, you have to write
the resulting message dictionaries to a file with json.dump() and add
a newline after each message.

.. code-block:: python

    import ais.stream
    import json
    import ais.compatibility.gpsd

    with open("infile.nmea") as inf:
        with open("outfile.gpsd") as outf:
            for msg in ais.stream.decode(f):
                gpsdmsg = ais.compatibility.gpsd.mangle(msg)
                json.dump(gpsdmsg, outf)
                outf.write("\n")

AIS Specification Documents
---------------------------

- ITU-1371, ITU-1371-{1,2,3,4]
- `e-Navigation <http://www.e-navigation.nl/asm>`_
- IMO Circ 236
- IMO Circ 289
- EU RIS

Developing
----------

The C++ code was switched to the Google style in November, 2012.
Indenting should be by 2 spaces.

http://google-styleguide.googlecode.com/svn/trunk/cpplint/

.. code-block:: console

    $ git clone https://github.com/schwehr/libais
    $ cd libais
    $ pip install -e .\[test\]
    $ py.test test ais --cov ais --cov-report term-missing
