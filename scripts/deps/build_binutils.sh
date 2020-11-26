#!/usr/bin/env bash

BASH_DEBUG="${BASH_DEBUG:=}"
([ "${BASH_DEBUG}" == "true" ] || [ "${BASH_DEBUG}" == "1" ]) && set -o xtrace
set -o errexit
set -o nounset
set -o pipefail
shopt -s dotglob
trap "exit" INT

echo "On Ubuntu you will need these dependencies:"
echo "sudo apt-get install sharutils texi2html"
echo ""

VERSION=2.35.1

NAME=binutils
SRC_DIR="${THIS_DIR}/../tmp"
SOURCE_DIR=${SRC_DIR}/${NAME}-${VERSION}
BUILD_DIR=${SOURCE_DIR}_build
INSTALL_DIR="${THIS_DIR}/../3rdparty/${NAME}"
ARCHIVE_FILE=${NAME}-${VERSION}.tar.xz
URL=http://ftp.gnu.org/gnu/${NAME}/${ARCHIVE_FILE}

mkdir -p ${SRC_DIR}
pushd ${SRC_DIR}

# Download package
if [ ! -f ${ARCHIVE_FILE} ]; then
  wget ${URL}
fi

# Create build directory, configure and build
mkdir -p ${BUILD_DIR}
pushd ${BUILD_DIR}

${SOURCE_DIR}/configure \
  --prefix="${INSTALL_DIR}" \
  --enable-gold \
  --enable-ld=default \
  --enable-plugins \
  --enable-shared \
  --disable-werror \
  --enable-64-bit-bfd \
  --with-system-zlib

make -j ${NUM_JOBS}
make install

popd
