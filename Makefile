# Make sure $JSSHELL points to your js shell binary in .profile or .bashrc
# Most targets use commands that need a js shell path specified
JSSHELL ?= $(error Specify a valid path to a js shell binary in ~/.profile: export JSSHELL=C:\path\js.exe or /path/js)

# Version number used in naming release files.
VERSION ?= $(error Specify a version for your release (e.g., VERSION=0.5))
FREQUENCY ?= 440

benchmark:
	${JSSHELL} -m -j -p -e 'var FREQUENCY=${FREQUENCY};' -f ./dist/dsp.js -f ./bench/bench.js -f ./bench/dft.js
	${JSSHELL} -m -j -p -e 'var FREQUENCY=${FREQUENCY};' -f ./dist/dsp.js -f ./bench/bench.js -f ./bench/fft.js
	${JSSHELL} -m -j -p -e 'var FREQUENCY=${FREQUENCY};' -f ./dist/dsp.js -f ./bench/bench.js -f ./bench/rfft.js
	${JSSHELL} -m -j -p -f dsp.js -f ./bench/bench.js -f ./bench/deinterleave.js

clean:
	rm -fr ./release
