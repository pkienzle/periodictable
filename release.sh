#!/bin/bash

# Usage: ./release.sh [sync|test|doc|notes|web|pypi]
#
# Build an official elements package release
#
# Releasing this package requires some setup in your local environment
#    nosetests and coverage package
#    sphinx and latex
#    mathjax with \AA extension
#    hudson server set up to build/test on windows/mac
#    reflectometry.org server key for updating docs
#    ~/.pypirc should be defined
#
# The patched MathJax (see below) needs to be symlinked into your
# doc/sphinx directory.

# The following is a minimal patch to MathJax to use the Angstrom symbol in TeX
# == MathJax/unpacked/jax/input/TeX/jax.js ==
#           // Ord symbols
#           S:            '00A7',
# +         AA:           '212B',
#           aleph:        ['2135',{mathvariant: MML.VARIANT.NORMAL}],
#           hbar:         '210F',

# Adapt the following to your own username/password for pypi, and get yourself
# added to the periodictable package.
# == ~/.pyirc ==
# [distutils]
# index-servers =
#     pypi
#   
# [pypi]
# username:...
# password:...

case "$1" in
sync)   step=0;;
test)   step=1;;
doc)    step=2;;
notes)  step=3;;
web)    step=4;;
pypi)   step=5;;
*)      step=0;;
esac


function ready() {
  stepname=$1; shift
  echo -n "$* [y/n] "
  read ans && test "$ans" != "y" \
    && echo Restart with ./release.sh $stepname && exit
}

if [ $step -le 0 ]; then
  echo === Version control status ===
  set -x
  git pull
  git status
  set +x
  ready sync Is the repository up to date?
fi

if [ $step -le 1 ]; then
  echo === Tests ===
  set -x
  python2.7 test.py -q --with-coverage
  python2.6 test.py -q
  set +x
  if false; then
    echo
    # Ask hudson build server if package is working on all platforms
    hudson_server="localhost:8080"
    elements_job="elements"
    url="http://$hudson_server/job/$elements_job/lastBuild"
    jsonurl="$url/api/json?depth=0"
    if curl --silent $jsonurl | grep -q SUCCESS ; then
        echo latest hudson build was successful
    else
        echo **** latest hudson build failed ... see $url
        firefox $url &
    fi
  fi
  ready test Are the tests okay?
fi

if [ $step -le 2 ]; then
  echo === Documentation ===
  (cd doc/sphinx && make clean html pdf)
  firefox doc/sphinx/_build/html/index.html >/dev/null 2>&1 & 
  evince doc/sphinx/_build/latex/PeriodicTable.pdf >/dev/null 2>&1 &
  ready doc Does the documentation build correctly?
fi

if [ $step -le 3 ]; then
  echo === Release notes ===
  rst2html README.rst > /tmp/README.html
  firefox /tmp/README.html >/dev/null 2>&1 &
  git log --format="%Cred%ad%Creset %s %Cred%an%Creset" --date=short
  version=$(grep __version__ periodictable/__init__.py | sed -e's/^.*= *//')
  echo *** Version is $version
  ready notes Are the release notes shown in the browser up to date?
fi

if [ $step -le 4 ]; then
  ready web: Push docs to the web?
  ssh reflectometry.org rm -r web/danse/docs/elements
  find doc/sphinx/_build/html | xargs chmod ug+rw
  find doc/sphinx/_build/html -type d | xargs chmod g+x
  rm -r doc/sphinx/_build/html/_static/MathJax
  (cd doc/sphinx/_build && scp -r html reflectometry.org:web/danse/docs/elements)
  ssh reflectometry.org ln -s /var/www/reflectometry/MathJax web/danse/docs/elements/_static
  ready web Documentation upload successful?
fi

if [ $step -le 5 ]; then
  ready Push package to pypi?
  git tag -a v$version -m "Release $version"
  git push --tags
  python setup.py sdist upload
  ready pypi Package upload successful?
fi

echo == All done! ==
