#!/usr/bin/bash

if [[ $0 == $BASH_SOURCE ]]; then
    echo "Please source this script."
    exit 1
fi

export MU2E_SEARCH_PATH="${MU2E_SEARCH_PATH}:${PWD}"
