#!/usr/bin/env bash

BASH_DEBUG="${BASH_DEBUG:=}"
([ "${BASH_DEBUG}" == "true" ] || [ "${BASH_DEBUG}" == "1" ]) && set -o xtrace
set -o errexit
set -o nounset
set -o pipefail
shopt -s dotglob
trap "exit" INT

echo "On Ubuntu you will need these dependencies:"
echo "sudo apt-get install build-essential git swig cmake libc++abi-dev libc++-dev libc++-helpers doxygen  python2.7-dev libedit-dev libncurses5-dev libxml2-dev libffi-dev ocaml ocaml-findlib"
echo ""

THIS_DIR=$(
  cd $(dirname "${BASH_SOURCE[0]}")
  pwd
)

if [[ $OSTYPE == "linux-gnu" ]]; then
  export NUM_JOBS="$(nproc)"
elif [[ $OSTYPE == "darwin"* ]]; then
  export NUM_JOBS="$(sysctl -n hw.ncpu)"
fi

NAME="llvm"
URL="https://github.com/llvm/llvm-project"
VERSION_DEFAULT="10.0.0"
GIT_TAG_PREFIX="llvmorg-"

VERSION=${1:-${VERSION_DEFAULT}}
VERSION_MAJOR="$(echo ${VERSION} | cut -d. -f1)"
NAME_SIMPLE=${NAME}-${VERSION_MAJOR}
BRANCH=${VERSION}
NAME=${NAME}-${BRANCH}

if [ "${BRANCH}" == "master" ] || [ -z "${BRANCH}" ]; then
  COMMIT_HASH=$(git ls-remote ${URL} | grep HEAD | cut -c-7)
else
  BRANCH="${GIT_TAG_PREFIX}${BRANCH}"
  COMMIT_HASH=$(git ls-remote ${URL} | grep -i "refs/tags/${BRANCH}^{}$" | cut -c-7)
fi

if [ ! -z "${COMMIT_HASH}" ]; then
  NAME=${NAME}-${COMMIT_HASH}
fi

SRC_DIR="${THIS_DIR}/../tmp"
SOURCE_DIR=${SRC_DIR}/${NAME}
BUILD_DIR=${SOURCE_DIR}_build
INSTALL_DIR="${THIS_DIR}/../3rdparty/${NAME_SIMPLE}"

if [ ! -d ${SOURCE_DIR} ]; then
  mkdir -p ${SRC_DIR}
  git clone --recursive --depth 1 -b ${BRANCH} ${URL} ${SOURCE_DIR}
fi

mkdir -p ${BUILD_DIR}
pushd ${BUILD_DIR}

unset CFLAGS CXXFLAGS CMAKE_C_FLAGS CMAKE_CXX_FLAGS
ADDITIONAL_INCLUDE_PATH="/usr/include/x86_64-linux-gnu"
export C_INCLUDE_PATH="${ADDITIONAL_INCLUDE_PATH}${C_INCLUDE_PATH:+:$C_INCLUDE_PATH}"
export CPLUS_INCLUDE_PATH="${ADDITIONAL_INCLUDE_PATH}${CPLUS_INCLUDE_PATH:+:$CPLUS_INCLUDE_PATH}"

cmake "${SOURCE_DIR}/llvm" \
  -DCMAKE_INSTALL_PREFIX="${INSTALL_DIR}" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_VERBOSE_MAKEFILE=0 \
  -DLLVM_TOOL_LLDB_BUILD=1 \
  -DLLVM_TOOL_LLD_BUILD=1 \
  -DLLVM_ENABLE_FFI=1 \
  -DLLVM_CCACHE_BUILD=1 \
  -DLLVM_ENABLE_LIBCXX=1 \
  -DLLVM_ENABLE_LLD=1 \
  -DLLVM_ENABLE_PROJECTS="clang;clang-tools-extra;libcxx;libcxxabi;libunwind;lld;lldb;openmp" \
  -DDEFAULT_CXX_STDLIB="libc++" \
  -DCMAKE_C_FLAGS="-I/usr/include/x86_64-linux-gnu" \
  -DCMAKE_CXX_FLAGS="-I/usr/include/x86_64-linux-gnu" \
  -DLLVM_BINUTILS_INCDIR="${THIS_DIR}/../3rdparty/binutils/include"

make -j ${NUM_JOBS}
make install

popd
