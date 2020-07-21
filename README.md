[![Build Status](https://travis-ci.org/bnlnpps/tpc-rs.svg?branch=master)](https://travis-ci.org/bnlnpps/tpc-rs)

**tpc-rs** provides a library to simulate the signal detected by a Time
Projection Chamber (TPC) and digitized by the readout electronics when charge
particles traverse its sensitive volume.

## Quick Start

Building the library from source is fairly straightforward. The only external
dependency is [ROOT](https://github.com/root-project/root) which is assumed to
be available on the system.

    git clone https://github.com/bnlnpps/tpc-rs.git
    cd tpc-rs && mkdir build && cd build
    cmake ../
    cmake --build . -- install


## Installation

Additional arguments to the `cmake` command can be provided to customize the
installation if the default ones do not meet your project/system requirements.
Here is a few commonly used CMake options:

    -DCMAKE_INSTALL_PREFIX=<path>

`<path>` can be absolute or relative. One can specify it to distinguish between
different build types and/or versions, e.g.
`/opt/local/tpc-rs-ABC-Debug-i686` or `${HOME}/tpc-rs-XYZ-Release-x86_64`.

    -DCMAKE_BUILD_TYPE=<build_type>

with `<build_type>` being one of `Debug`, `RelWithDebInfo`, or `Release`.

     -DBUILD_TESTING=OFF

skips building of the tests.


## How to Use

To use an installed tpc-rs library in another CMake project just include the
following in your CMakeLists.txt:

    find_package(tpcrs [major.minor] [EXACT] [REQUIRED])
    target_link_libraries(<your_target> PUBLIC tpcrs)


## Basic Test

Once the source code is built one can execute basic tests to verify the
performance of the library. A reference file is provided with a pre-generated
input and the corresponding output from the main `tpc-rs` conversion routine.

    cd tpc-rs/build
    python ../utils/gdrive.py 1XESdyqg6kXhu5gPrgS3drq3UIBu0f1l5 | tar -xz
    ../tests/run_tests.sh
