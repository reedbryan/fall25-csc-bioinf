# Bioinformatics with Codon

This project implements a De Bruijn Graph (DBG) assembler and N50 calculator using Codon.

## Setup and Installation

1. Install Codon:
```bash
/bin/bash -c "$(curl -fsSL https://exaloop.io/install.sh)"
```

2. Install the Seq plugin for bioinformatics:
```bash
curl -L https://github.com/exaloop/seq/releases/download/v0.11.5/seq-$(uname -s | tr '[:upper:]' '[:lower:]')-$(uname -m).tar.gz | tar zxvf - -C ${HOME}/.codon/lib/codon/plugins
```

3. Test your installation:
```bash
echo -ne "import bio\nprint('hello', s'ACGT')\n" | ~/.codon/bin/codon run -plugin seq -
```

## Running the Code

To run the De Bruijn Graph assembler and calculate N50:

```bash
~/.codon/bin/codon run -plugin seq main.codon data1
```

Replace `data1` with the name of your data directory.

## VS Code Integration

For proper syntax highlighting and code completion:

1. VS Code should treat `.codon` files as Python files (configured in `.vscode/settings.json`)
2. When running Codon files, always use the `-plugin seq` flag to enable bioinformatics functionality

## Troubleshooting

If you're seeing "Import could not be resolved" or similar errors in VS Code, this is normal as the IDE doesn't have native support for Codon's bio libraries. The code will still run correctly when executed with the Codon interpreter and seq plugin.
