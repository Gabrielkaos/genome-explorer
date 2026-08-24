export const BASE_COLORS = ["#68c98f", "#5ec8d8", "#e8c15a", "#ef7fa3"];
export const GAP_COLOR = "#39424e";
export const BASE_LETTERS = ["A", "C", "G", "T", "-"];

export function baseColor(code) {
  return code > 3 ? GAP_COLOR : BASE_COLORS[code];
}
