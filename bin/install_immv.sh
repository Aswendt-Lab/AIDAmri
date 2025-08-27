#!/bin/bash 

immv_wrapper_code=$(cat <<'EOF'
# replace the first line with your fslpython path
# -*- coding: utf-8 -*-
import re
import sys

from fsl.scripts.immv import main

if __name__ == '__main__':
    sys.argv[0] = re.sub(r'(-script\.pyw?|\.exe)?$', '', sys.argv[0])
    sys.exit(main())
EOF
)

# Create the immv wrapper called by bet4animal
if [ -f "${FSLDIR}/bin/immv" ];
  then
    echo "immv already exists"
else
    echo "creating immv file interpreter line"       
    echo "#!"`ls ${FSLDIR}/fslpython/bin/python` > ${FSLDIR}/bin/immv
    # Add the immv function code
    echo "${immv_wrapper_code}" >> ${FSLDIR}/bin/immv
    # Make the immv wrapper executable
    chmod +x ${FSLDIR}/bin/immv
fi
