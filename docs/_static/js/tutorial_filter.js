document.addEventListener("DOMContentLoaded", function() {
    const filterButtons = document.querySelectorAll(".filter-btn");
    const cards = document.querySelectorAll(".tutorial-card");

    let activeFilters = new Set();

    filterButtons.forEach(button => {
        button.addEventListener("click", () => {
            const filter = button.dataset.filter;

            if (filter === "all") {
                activeFilters.clear();
            } else {
                if (activeFilters.has(filter)) {
                    activeFilters.delete(filter);
                } else {
                    activeFilters.add(filter);
                }
            }

            // Toggle "active" class for styling
            button.classList.toggle("active");

            // Display cards based on active filters
            cards.forEach(card => {
                const cardTags = card.dataset.tags.split(", ");
                if (activeFilters.size === 0 || [...activeFilters].every(f => cardTags.includes(f))) {
                    card.style.display = "";
                } else {
                    card.style.display = "none";
                }
            });
        });
    });
});
