/* global describe it */
var assert = require("assert");
var dsp = require("..");

describe("dsp.js library", function () {
  it("Export all modules", function () {
    var keys = Object.keys(dsp);
    assert.equal(16, keys.length, "must export 16 items");
    assert(typeof dsp["DSP"] === "object", "DSP must be an object");
    keys.forEach(function (key) {
      if (key !== "DSP") {
        assert(typeof dsp[key] === "function", key + " must be a function.");
      }
    });
  });
});
