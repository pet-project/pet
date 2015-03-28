
# PET (Paral's Plasma Physics Toolkit) #

PET project for plasma simulations

## Prerequisites ##

* [Cython](http://www.cython.org) (>= 0.19)
* sphinx (>= 1.3)
* sphinx-bootstrap-theme (>= 0.4):
  [GitHub](http://github.com/mcmtroffaes/sphinxcontrib-bibtex/),
  [Documentation](http://sphinxcontrib-bibtex.readthedocs.org/en/latest/)
* sphinxcontrib-bibtex (>= 0.3):
  [GitHub](http://github.com/mcmtroffaes/sphinxcontrib-bibtex/),
  [Documentation](http://sphinxcontrib-bibtex.readthedocs.org/)
  
## Download ##

For read only access:

```bash
cd
git clone https://github.com/pet-project/pet.git
```

For development, you must set up your GitHub account with SSH keys and checkout
`PET` as follows:

```bash
cd
git clone git@github.com:pet-project/pet.git
```

Setting up SSH keys and using native GIT transport protocol is better in that
you don't have to provide credentials all the time.

## Building the library ##

You can build the Python bindings for PET by calling *make*

## Development Environment ##

If you would like to set up a development environment, the easiest way is to
source a shell script which sets up all variables for you. This script requires
only `PET_ROOT` variable to be set beforehand. Here is an example how to do
that assuming that `bash` is your default shell. Modify your `$HOME/.bashrc`
file by adding following lines (modify as needed):

```bash
export PET_HOME="$HOME/pet"
source "$PET_HOME/bin/env.sh"
```

## Pull Request ##

Please submit a pull reguest through GitHub interface.

If you are submitting a major change which you would like to mention in change
log file, include a new file with changes in the same format as `CHANGELOG.md`
and place the file in `doc/incoming/name_of_feature.md`. This is particularly
useful when you are submitting a change which is not backward compatible or
some major enhancement. Once we are ready for release we will move the changes
into `CHANGELOG.md` file.

## Contact Us ##

Please contact us at any of these communication channels:
* [PET-usr](https://groups.google.com/forum/#!forum/pet-usr) Google group for
  users.
* [PET-dev](https://groups.google.com/forum/#!forum/pet-dev) Google group for
  developers.
