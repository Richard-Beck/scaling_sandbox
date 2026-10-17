# Development output organization

Group each experiment and all its exports in its topic folder; keep only this index, these instructions and experiment folders at dev/. Use the existing supporting/{vectors,data,metadata,viewers} structure inside each experiment.

After regeneration or adding a new export, run code/mechanisms/organize_length_dev.py from the manuscript directory (or its repository-relative path). Update the prefix routing there for a new experiment. Preserve raw simulation evidence and frozen campaign sources. Verify local links and artifact hashes after moves.
