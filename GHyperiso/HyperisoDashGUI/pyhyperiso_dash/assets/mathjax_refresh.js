(function () {
  function collectTargets() {
    var selectors = [
      '.dash-table-container',
      '.js-plotly-plot',
      '.MathJax-auto-typeset',
      '.mathjax-auto-typeset'
    ];
    var nodes = [];
    selectors.forEach(function (sel) {
      document.querySelectorAll(sel).forEach(function (el) {
        if (!el.closest('.Select-menu-outer') && !el.closest('.Select__menu')) {
          nodes.push(el);
        }
      });
    });
    return nodes;
  }

  function typeset() {
    if (window.MathJax && window.MathJax.typesetPromise) {
      window.MathJax.typesetPromise(collectTargets())
        .catch(function (err) { console.warn('MathJax typeset failed:', err); });
    }
  }

  var timer = null;
  function scheduleTypeset() {
    if (timer) {
      window.clearTimeout(timer);
    }
    timer = window.setTimeout(typeset, 160);
  }

  document.addEventListener('DOMContentLoaded', function () {
    scheduleTypeset();
    var target = document.body;
    if (!target || !window.MutationObserver) {
      return;
    }
    var observer = new MutationObserver(function (mutations) {
      for (var i = 0; i < mutations.length; i++) {
        if (mutations[i].addedNodes && mutations[i].addedNodes.length) {
          scheduleTypeset();
          return;
        }
      }
    });
    observer.observe(target, { childList: true, subtree: true, characterData: true });
  });
})();
