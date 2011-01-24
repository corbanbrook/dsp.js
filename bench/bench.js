function benchmark(func, loopCount) {
  loopCount = loopCount || 10000;

  var start = Date.now();

  for (var i = 0; i < loopCount; i++) { 
    func();
  }

  var end = Date.now();
  return end - start;
}
