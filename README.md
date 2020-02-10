**tpcrs** is a library to simulate response from a Time Projection Chamber (TPC)
detector when it collects and digitizes the signal produced by charge particles
traversing its volume.

## Quick Start

Building the library from source is straightforward. The only external
dependency is [ROOT](https://github.com/root-project/root) which is assumed to
be installed and available.

    git clone https://github.com/plexoos/tpc-rs.git
    cd tpc-rs && mkdir build && cd build
    cmake ../
    make

## Basic Test

Once the source code is built one can execute a basic test to verify the
performance of the library. A reference file is provided with a pre-generated
input and the corresponding output from the main `tpcrs` conversion routine.

    cd tpc-rs/build
    python ../utils/gdrive.py 13GZLh1XClfJb-sQkhTxG9Rfe7owOYcZE > geant_event.root
    export STAR=/path/to/tpc-rs/
    ./tests/test_tpcrs
    diff test_tpcrs_inp.log test_tpcrs_out.log
