# RoundingBack
Repository for RoundingBack, a native backbone extractor for Pseudo-Boolean Optimization based on RoundingSat.

To extract the backbone for an instance you first need to run in the root directory of roundingback:

    cd build
    cmake -DCMAKE_BUILD_TYPE=Release ..
    make

Once it is built, you need to run:

    ./extract_backbone.sh <pbo_instance_path>

You can also pass a timeout and/or an optimal value for an optimization instance if you already know it:

    ./extract_backbone.sh <pbo_instance_path> --opt <optimum_value> --timeout <time_in_seconds>
