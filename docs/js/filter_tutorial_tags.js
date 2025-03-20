// filtering process inspired by PyTorch Lightning sphinx theme

window.filterTags = {
    bind: function() {
      const ALL_TAG_CLASS = "all-tag-selected";
      const SELECTED_CLASS = "selected";
      const FILTER_BTN_SELECTOR = ".filter-btn";
      const ALL_TAG_SELECTOR = `[data-tag='all']`;
  
      const options = {
        valueNames: [{ data: ["tags"] }],
        page: 10,
        pagination: true
      };
      
      const tutorialList = new List("tutorial-cards", options);
  
      // Check if a tutorial's tags match the selected tags
      const tagsMatch = (cardTags, selectedTags) => 
        selectedTags.every(tag => cardTags.includes(tag));
  
      // Filter the list based on selected tags
      const filterTutorials = () => {
        const selectedTags = Array.from(document.querySelectorAll(`.${SELECTED_CLASS}`))
                                  .map(el => el.dataset.tag);
  
        tutorialList.filter(item => {
          const cardTags = (item.values().tags ?? "").split(",");
          return selectedTags.length === 0 || tagsMatch(cardTags, selectedTags);
        });
      };
  
      // Event listener for filter button clicks
      document.querySelectorAll(FILTER_BTN_SELECTOR).forEach(button => {
        button.addEventListener("click", () => {
          const tag = button.dataset.tag;
  
          if (tag === "all") {
            button.classList.add(ALL_TAG_CLASS);
            document.querySelectorAll(`.${SELECTED_CLASS}`).forEach(el => el.classList.remove(SELECTED_CLASS));
          } else {
            button.classList.toggle(SELECTED_CLASS);
            document.querySelector(ALL_TAG_SELECTOR).classList.remove(ALL_TAG_CLASS);
          }
  
          if (!document.querySelector(`.${SELECTED_CLASS}`)) {
            document.querySelector(ALL_TAG_SELECTOR).classList.add(ALL_TAG_CLASS);
          }
  
          filterTutorials();
        });
      });
    }
  };
  