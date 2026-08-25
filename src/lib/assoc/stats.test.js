import { describe, it, expect } from 'vitest';
import { fisherExactTwoSided, oddsRatioCI } from './stats.js';

describe('stats.js', () => {
  it('fisherExactTwoSided should handle 2x2 contingency tables', () => {
    const p = fisherExactTwoSided(2, 0, 0, 2);
    expect(p).toBeCloseTo(0.3333333, 5); 
  });

  it('oddsRatioCI should compute OR and 95% CI with Haldane-Anscombe correction', () => {
    const res = oddsRatioCI(10, 2, 3, 15);
    expect(res.or).toBeGreaterThan(1);
    expect(res.lo).toBeLessThan(res.or);
    expect(res.hi).toBeGreaterThan(res.or);
  });
});
