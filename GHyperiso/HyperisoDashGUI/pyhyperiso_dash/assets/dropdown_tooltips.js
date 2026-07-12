(function () {
  "use strict";

  const SELECTOR = [
    ".Select-value-label",
    ".Select-option",
    ".VirtualizedSelectOption",
    ".Select__single-value",
    ".Select__multi-value__label",
    ".Select__option"
  ].join(",");

  function bestTitle(node) {
    const explicit = node.querySelector && node.querySelector("[title]");
    if (explicit && explicit.getAttribute("title")) {
      return explicit.getAttribute("title").trim();
    }
    return (node.textContent || "").replace(/\s+/g, " ").trim();
  }

  function decorate(root) {
    const nodes = [];
    if (root && root.nodeType === 1 && root.matches && root.matches(SELECTOR)) {
      nodes.push(root);
    }
    if (root && root.querySelectorAll) {
      nodes.push(...root.querySelectorAll(SELECTOR));
    }
    nodes.forEach((node) => {
      const title = bestTitle(node);
      if (title) {
        node.setAttribute("title", title);
        node.setAttribute("aria-label", title);
      }
    });
  }

  function start() {
    decorate(document.body);
    const observer = new MutationObserver((mutations) => {
      mutations.forEach((mutation) => {
        mutation.addedNodes.forEach(decorate);
      });
    });
    observer.observe(document.body, { childList: true, subtree: true });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", start, { once: true });
  } else {
    start();
  }
})();
