/**
 * Returns the hyperbolic sine of the number
 *
 * @meta version: 1004.2314
 * @meta discuss at: http://phpjs.org/functions/sinh
 * @meta original by: Onno Marsman
 *
 * @param {Number} num
 * @example
 * sinh(-0.9834330348825909); // => -1.1497971402636502
 */
function sinh (arg) {
  return (Math.exp(arg) - Math.exp(-arg))/2;
}

module.exports = sinh;
