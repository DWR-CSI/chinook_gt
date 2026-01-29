#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process TEST_MOUNT {
    tag "Test Azure File Share"
    container 'quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:bd996097b6dd518cf788ddd6c586fb23d039cb9c-0'

    script:
    """
    echo "Listing contents of /mnt/batch/tasks/fsmounts/seqres:"
    ls -la /mnt/batch/tasks/fsmounts/seqres || echo "Failed to list directory"
    
    echo "Listing contents of /mnt/seqres (checking alternate path):"
    ls -la /mnt/seqres || echo "Failed to list directory"
    
    echo "Mount verification complete."
    """
}

workflow {
    TEST_MOUNT()
}
