#!/bin/bash

echo "Setting up...\n"

# If snp-pileup is not in user PATH, edit this
SNP_PILEUP="snp-pileup"

sed 's/snp-pileup/${SNP_PILEUP}/' snp-pileup-wrapper.R

echo "Done"