
---

# RNAProbesSource

**RNAProbesSource** is the source code for the RNAProbes suite. This repository is a subtree of RNAProbes and is intended to be treated as read-only.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)  
  - [Using pip](#using-pip)  
  - [Using Poetry](#Using-Poetry-recommended-for-linux)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

---

## Installation

### Clone the repository
```bash
git clone https://github.com/Catrina-Lab/RNAProbesSource.git
```

### Using pip

#### Windows, macOS, or linux with venv with `pip`:
```bash
cd RNAProbesSource
pip install -r requirements.txt
````

### Using Poetry (recommended for linux)

Ensure you have [Poetry](https://python-poetry.org/docs/#installation) installed (make sure to follow the instructions there to install poetry >= 2.0).

#### Clone the repo and install dependencies:

```bash
poetry install --project RNAProbesSource #replace RNAProbesSource with the directory name
```

#### Running the package inside Poetry’s virtual environment:

```bash
poetry run --project RNAProbesSource python -m RNAProbesSource
```

---

## Usage

Once installed, you can run the main program:

```bash
cd .. #make sure you're in the directory above RNAProbesSource
python -m RNAProbesSource
```

Or, if using Poetry:

```bash
cd .. #make sure you're in the directory above RNAProbesSource
poetry run --project RNAProbesSource python -m RNAProbesSource
```

---

## Contributing

We welcome any contributions. However, this repository is a readonly subtree. To contribute:
1. Fork the main [RNAProbes repository](https://github.com/Catrina-Lab/RNAProbes)
2. Clone your fork and make desired changes.
4. Submit a pull request against the main RNAProbes repo. If accepted, it will automatically be changed here also.

Pull requests should include:

* Clear description of changes
* Updated tests (if applicable)
* Documentation updates (as needed)
---

## License

GNU GPL 3.0

---


