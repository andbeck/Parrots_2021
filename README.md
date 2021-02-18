# Parrots_2021

## Parrot Demographic Modelling

To get started with this repository, clone it to a suitable location from the command line: 

```
cd <target_directory>
git clone https://github.com/andbeck/Parrots_2021.git
```

## Code structure

### `data/`

- Contains 2006-2014 spreadsheet focused on Bonaire long-term collection
- Note that [issue #1](https://github.com/andbeck/Parrots_2021/issues/1) from @tamora contains further data sources from YSA and other species - such as Lilac Crowned Amazon - this needs to be reviewed
- Need to generate table like in appendices in Wisdom - a table with low uncertainty

### `doc/`

- Contains starter for Markdown Files for reporting

### `fig/`

- Contains picture for markdown reporting highlighting the life cycle currently used
- Could also be used to save figures from analyses

### `R/`

- `main.R` Initial outline for code to analyse multiple matrices - needs LSA
- `make_projection_matrix.R` Code to construct matrices
- `YSA_life_history_data.R` Code to construct mean and se from YSA core data from Bonaire 2004-2016
