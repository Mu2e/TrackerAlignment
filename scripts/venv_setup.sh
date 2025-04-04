#!/usr/bin/bash

if [[ $0 == $BASH_SOURCE ]]; then
    echo "Please source this script."
    exit 1
fi

if [ ! -f "${TRKALIGN_BASE}/.no-venv" ]; then
    if [ ! -d "${TRKALIGN_BASE}/.venv" ]; then
        echo "Creating virtual environment at ${TRKALIGN_BASE}/.venv"
        python -m venv ${TRKALIGN_BASE}/.venv 
        source ${TRKALIGN_BASE}/.venv/bin/activate

        python -m pip install -r ${TRKALIGN_BASE}/scripts/requirements.txt
    else
        echo "Sourcing virtual environment at ${TRKALIGN_BASE}/.venv"
        source ${TRKALIGN_BASE}/.venv/bin/activate
    fi
else
    echo "Skipped virtual environment setup."
fi
# set up some convenience commands 

alias aligntrack_display='python ${TRKALIGN_SCRIPTS_DIR}/aligntrack_display.py ' 

