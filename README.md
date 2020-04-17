**tpc-rs** provides a library to simulate the signal detected by a Time
Projection Chamber (TPC) and digitized by the readout electronics when charge
particles traverse its sensitive volume.

## Quick Start

Building the library from source is straightforward. The only external
dependency is [ROOT](https://github.com/root-project/root) which is assumed to
be installed on the system.

    git clone https://github.com/plexoos/tpc-rs.git
    cd tpc-rs && mkdir build && cd build
    cmake ../
    make

Note that ROOT version less than 6.18 is required due to dependecy on the Table
module removed from the more recent releases. We are working on removing this
limitation for `tpc-rs`.

## Basic Test

Once the source code is built one can execute basic tests to verify the
performance of the library. A reference file is provided with a pre-generated
input and the corresponding output from the main `tpc-rs` conversion routine.

    cd tpc-rs/build
    python ../utils/gdrive.py 1XESdyqg6kXhu5gPrgS3drq3UIBu0f1l5 | tar -xz
    ../tests/run_tests.sh
