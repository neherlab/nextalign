<h1 id="nextclade" align="center">
Nextalign
</h1>

<h4 id="nextclade" align="center">
Sequence alignment
</h1>


## Developer's guide

### Quick start

 1. Install and configure required dependencies
 
    - cmake >= 3.10
 
    - conan package manager (or Ubuntu Linux there is an installation script included in `./tools/install-conan`)
 
    - nodemon (requires Node.js)
    
    ```bash
    npm install --global nodemon
    ```


 2. Clone and run

    ```bash
    git clone --recursive https://github.com/neherlab/nextalign
    cd nextalign
    make dev

    ```
    (note the `--recursive` flag for git - this repositoy contains git submodules)

    This will:
    
     - install or update conan packages
     - run cmake and generate makefiles
     - build the project and tests
     - run static analysis on source files
     - run tests
     - watch source files and rebuild on changes
     
     If you don't want to install Node.js and nodemon, or don't want the automatic rebuild, you can use `make dev-nowatch` instead of `make dev`.


### Tests

Test are run as a part of the main development script (`make dev`). We are using [Google Test](https://github.com/google/googletest/blob/master/googletest/docs/primer.md)
See [Google Test documentation](https://github.com/google/googletest/blob/master/googletest/docs/primer.md) and [Google Mock documentation](https://github.com/google/googletest/blob/master/googlemock/README.md) for more details.

## License

<a target="_blank" rel="noopener noreferrer" href="LICENSE" alt="License file">MIT License</a>
