#!/bin/sh

aliroot -b << EOF
gSystem->Load("libMUONrec");
.x PatchHVPbPb2018.C+
.q
EOF
