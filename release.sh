#!/bin/bash

# change to the dir of the script
cd $( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

unittest() {
    python -m unittest test.test_suite_all
}

increase_patch() {
    git describe --tags --abbrev=0 |
    awk -F "." '
        {
            printf("%s.%d.%d", $1, $2, $3 + 1)
        }
    '
}

increase_minor() {
    git describe --tags --abbrev=0 |
    awk -F "." '
        {
            printf("%s.%d.0", $1, $2 + 1)
        }
    '
}

increase_major() {
    git describe --tags --abbrev=0 |
    awk -F "." '
        {
            printf("v%d.%d.%d", substr($1, 2) + 1, $2, $3)
        }
    '
}

release_github() {
    git commit -am "release ${version}"
    git push
    gh release create "${version}" --notes "release ${version}"
    git pull
}

update_bioconda_meta() {
    local tarball=${url##*/}
    local num_version=$(sed -r 's/^v//' <<<${version})
    local escaped_version=$(sed -r 's/\./\\\./g' <<<${num_version})
    local url_with_dynamic_version=$(sed -r 's/'${escaped_version}'/{{ version }}/' <<<${url})
    if ! [ -f "${version}.tar.gz" ]
    then
        wget $url
    fi
    local sha256=$(sha256sum ${tarball} | cut -d' ' -f1)

    sed -r \
        -e '/^\{% set version = ".*" %\}$/ s/"(.*)"/"'"${num_version}"'"/' \
        -e '/^  url: / s|^(  url: )(.*)$|\1'"${url_with_dynamic_version}"'|' \
        -e '/^  sha256: / s|^(  sha256: )(.*)$|\1'"${sha256}"'|' \
        "deploy/bioconda/meta.yaml.template"
}

test_bioconda() {
    mkdir -p ${BIOCONDA_RECIPES}/recipes/${pkg}
    # Update files
    cp LICENSE.md ${BIOCONDA_RECIPES}/recipes/${pkg}/
    update_bioconda_meta \
        > ${BIOCONDA_RECIPES}/recipes/${pkg}/meta.yaml
    pushd ${BIOCONDA_RECIPES}
    conda build recipes/${pkg}
    popd
}

release_bioconda() {
    pushd ${BIOCONDA_RECIPES}

    git checkout --force master
    # Delete local branch
    git branch -D "update_${pkg}"
    # Delete branch in your fork via the remote named "origin"
    git push origin -d "update_${pkg}"

    # exit on error
    set -e

    # Make sure our master is up to date with Bioconda
    git pull upstream master
    git push origin master

    # Create and checkout a new branch
    git checkout -b "update_${pkg}"

    popd

    # Update files
    cp LICENSE.md ${BIOCONDA_RECIPES}/recipes/${pkg}/
    update_bioconda_meta \
        > ${BIOCONDA_RECIPES}/recipes/${pkg}/meta.yaml

    pushd ${BIOCONDA_RECIPES}

    git commit -am "Update ${pkg}"

    git push --set-upstream origin "update_${pkg}"

    gh pr create --repo bioconda/bioconda-recipes --fill --template PULL_REQUEST_TEMPLATE.md
    popd
}

# action: github (create new release on github), test (test the latest github release for bioconda), deploy (deploy the latest github release to bioconda)
action=$1

if [[ "${action}" =~ "github" ]]
then
    version="$(increase_patch)"
else
    version="$(git describe --tags --abbrev=0)"
fi
pkg="rearr"
url="https://github.com/ljw20180420/${pkg}/archive/refs/tags/${version}.tar.gz"

unittest
if [[ "${action}" =~ "github" ]]
then
    release_github
fi
if [[ "${action}" =~ "test" ]]
then
    test_bioconda
fi
if [[ "${action}" =~ "deploy" ]]
then
    release_bioconda
fi
