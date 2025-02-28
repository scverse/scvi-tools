function highlightCategory(category) {
    // Remove highlight from all rows first
    document.querySelectorAll('.tutorial-row').forEach(row => {
        row.classList.remove('highlight');
    });
    // Highlight only the relevant category
    document.querySelectorAll('.' + category).forEach(row => {
        row.classList.add('highlight');
    });
}
