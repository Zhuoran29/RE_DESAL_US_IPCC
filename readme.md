# Renewable Energy Powered Desalination in the US under climate change scenarios

## Directory of this repositoy
- Wind_resource: Includes API to download wind hourly data from WindToolkit
- VAGMD_batch  : Includes python script and data for the batch-VAGMD model


## This program relies on a certain version of WaterTAP-REFLO packages, please install the package in order to run the models:

**WaterTAP-REFLO** supports Python versions 3.8 through 3.10.

### Prerequisites

- The conda package and environment manager, for example by using the [Miniconda installer](https://docs.conda.io/en/latest/miniconda.html#miniconda) following the steps appropriate for your operating system

### Installation

To install **WaterTAP-REFLO**, run:

```sh
# Direct to the folder 'watertap_reflo'

cd watertap_reflo
conda create --yes --name watertap-reflo-dev python=3.10 && conda activate watertap-reflo-dev
pip install -r requirements-dev.txt
```

### Running tests

```sh
conda activate watertap-reflo-dev
pytest --pyargs watertap_contrib.reflo
```


