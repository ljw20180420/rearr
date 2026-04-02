#!/bin/bash

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

version=$(increase_patch)
git commit -am "release ${version}"
git push
gh release create "${version}" --notes "release ${version}"
git pull
