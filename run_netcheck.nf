#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process NETCHECK {
    container 'python:3.11-slim'
    cpus 1
    memory '2 GB'
    publishDir "${launchDir}/results/netcheck", mode: 'copy'

    output:
    path 'netcheck.txt'

    script:
    """
    set +e
    {
      echo "--- python version ---"
      python3 --version
      echo "--- urllib github ---"
      python3 -c "
import urllib.request, traceback
try:
    with urllib.request.urlopen('https://api.github.com', timeout=15) as r:
        print('STATUS', r.status)
except Exception as e:
    traceback.print_exc()
"
      echo "--- DNS ---"
      python3 -c "import socket; print(socket.gethostbyname('api.github.com'))" 2>&1
      echo "--- env ---"
      env | grep -i github || echo "no GITHUB env vars set"
    } > netcheck.txt 2>&1
    cat netcheck.txt
    """
}

workflow {
    NETCHECK()
}
