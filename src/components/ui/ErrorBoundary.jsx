import React from "react";
import { C, FONT } from "../../theme.js";
import { Panel, Eyebrow } from "./Primitives.jsx";

export class ErrorBoundary extends React.Component {
  constructor(props) {
    super(props);
    this.state = { hasError: false, error: null };
  }

  static getDerivedStateFromError(error) {
    return { hasError: true, error };
  }

  componentDidCatch(error, errorInfo) {
    console.error("ErrorBoundary caught an error:", error, errorInfo);
  }

  render() {
    if (this.state.hasError) {
      return (
        <Panel style={{ padding: 20, borderLeft: `2px solid ${C.bad}` }}>
          <Eyebrow color={C.bad}>Runtime Error</Eyebrow>
          <div style={{ fontFamily: FONT, color: C.text, fontSize: 13, marginBottom: 12 }}>
            Something went wrong while rendering this section.
          </div>
          <div style={{
            background: "#05070a", padding: 12, borderRadius: 2,
            border: `1px solid ${C.border}`, color: C.bad,
            fontFamily: FONT, fontSize: 11, whiteSpace: "pre-wrap", overflowX: "auto"
          }}>
            {this.state.error && this.state.error.toString()}
          </div>
          <button onClick={() => this.setState({ hasError: false, error: null })}
            style={{
              all: "unset", cursor: "pointer", marginTop: 12, fontSize: 12, color: C.textDim,
              border: `1px solid ${C.border}`, padding: "4px 12px", borderRadius: 2, fontFamily: FONT
            }}>
            Try Again
          </button>
        </Panel>
      );
    }
    return this.props.children;
  }
}
