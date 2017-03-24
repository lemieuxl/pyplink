#!/bin/bash -x

# Getting pyplink's version to build
pyplink_version=$1
if [ -z $pyplink_version ]; then
    (>&2 echo "usage: $0 VERSION PYTHON_VERSION")
    exit 1
fi

# Getting the version of python we want to build for
python_version=$2
if [ -z $python_version ]; then
    (>&2 echo "usage: $0 VERSION PYTHON_VERSION")
    exit 1
fi

# Creating a directory for the build module
mkdir -p conda_dist

# Creating a directory for the skeleton
mkdir -p skeleton
pushd skeleton

# Creating the skeleton
conda skeleton pypi pyplink --version $pyplink_version

# Checking that fetching pyplink was successful
if [ $? -ne 0 ]; then
    (>&2 echo "Error when creating skeleton for pyplink version $pyplink_version")
    exit 1
fi

# The different python versions build
if [ "$python_version" == "2" ]; then
    python_versions="2.7"
elif [ "$python_version" == "3" ]; then
    python_versions="3.3 3.4 3.5 3.6"
else
    (>&2 echo "$python_version: invalid Python version")
    exit 1
fi

# The different build platforms
platforms="all"

# Building
for python_build_version in $python_versions; do
    # Building
    conda build --python $python_build_version pyplink &> log.txt

    # Checking the build was completed
    if [ $? -ne 0 ]; then
        cat log.txt
        (>&2 echo "Error when building pyplink $pyplink_version (python $python_build_version)")
        exit 1
    fi

    # Fetching the file name of the build
    filename=$(grep -oP "anaconda upload \K(\S+)$" log.txt)

    # Checking the file exists
    if [ -z $filename ]||[ ! -e $filename ]; then
        echo "Problem fetching file $filename" 1>&2
        exit 1
    fi

    # Converting to the different platforms
    for platform in $platforms; do
        conda convert -p $platform $filename -o ../conda_dist

        # Checking the conversion was completed
        if [ $? -ne 0 ]; then
            (>&2 echo "Problem converting pyplink $pyplink_version (python $python_build_version) to $platform")
            exit 1
        fi

    done

    # Purging
    conda build purge

done

popd
rm -rf skeleton

# Indexing
pushd conda_dist
conda index *
popd
