#!/bin/bash

function log_error () {
    echo -e "\e[31m -- ERROR! $1 \e[0m"
    exit 1
}

function log_success () {
    echo -e "\e[32m -- $1 \e[0m"
}

function log_message () {
    echo -e " -- $1"
}

if [ -z "$UCESB_DIR" ]; then
    log_error "UCESB_DIR environment variable is not set"
fi

if [ -z "$UCESB_BASE_DIR" ]; then
    log_error "UCESB_BASE_DIR environment variable is not set"
fi

if [ -z "$UPEXPS_DIR" ]; then
    log_error "UPEXPS_DIR environment variable is not set"
fi

if [ -z "$FAIRROOTPATH" ]; then
    log_error "FAIRROOTPATH environment variable is not set."
fi

if [ -z "$SIMPATH" ]; then
    log_error "SIMPATH environment variable is not set."
fi

if [ -z "$LMD_PATH" ]; then
    log_error "LMD_PATH environment variable is not set."
fi

if [ -z "$OUTROOTPATH" ]; then
    log_error "OUTROOTPATH environment variable is not set."
fi

if [ -z "$PARAMDIR" ]; then
    log_error "PARAMDIR environment variable is not set."
fi

if [ -z "$VMCWORKDIR" ]; then
    log_error "VMCWORKDIR environment variable is not set."
fi

log_success "All environment variables are set!"

