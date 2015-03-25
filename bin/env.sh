
#
# Modify your `$HOME/.bashrc` file by adding following lines (modify as needed)
#
# export PET_HOME="$HOME/pet"
# source "$PET_HOME/bin/env.sh"
#

if [ ! -d "$PET_HOME" ]; then
    echo ""
    echo "****  PET_HOME Variable not set! ****"
    echo ""
    return
fi

export PYTHONPATH="$PYTHONPATH:$PET_HOME"
