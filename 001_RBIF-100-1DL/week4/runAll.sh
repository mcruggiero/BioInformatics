#!/bin/bash

echo "removing old reports"
rm report.txt
./copyExomes.sh
./createCrisprReady.sh
./identifyCrisprSite.sh
./editGenome.sh
echo "completed"
