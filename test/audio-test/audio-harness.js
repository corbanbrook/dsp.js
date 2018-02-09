load('samples.js');

var startTime,
    totalTime;

function calcTime() {
  totalTime = (new Date()).getTime() - startTime;
}

function printResults(iterations) {
  print('Total Time: ' + totalTime + 'ms for ' + iterations + ' iterations, ' +
        (totalTime / iterations) + 'ms per iteration.');
}

function runTest(test, iterations) {
 startTime = new Date().getTime();
 for (var i = 0; i < iterations; i++) {
   test();
 }
 calcTime();
 printResults(iterations)
}
