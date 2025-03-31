// Create a set of all unique tags that appear in the cards
function getUniqueTags() {
    let tagSet = new Set();

    $(".card-container").each(function() {
        let tags = $(this).data('tags').split(",")
        .map(tag => tag.trim()) // Get rid of whitespace
        .filter(tag => tag !== ""); // Make sure there are no empty tags

        tags.forEach(tag => tagSet.add(tag)); // Add tag to the set
    });

    let uniqueTags = Array.from(tagSet);
    uniqueTags.sort();
    return uniqueTags;
}

// Get all of the tutorial group names
function getUniqueGroups() {
    let groupSet = new Set();

    $(".card-container").each(function() {
        let group = $(this).data('model-group-name')

        if (group && group.trim() !== "") {
            groupSet.add(group.trim());
        }

        groupSet.add(group)
    });

    let uniqueGroups = Array.from(groupSet);
    uniqueGroups.sort();
    return uniqueGroups;
}

// Create the tag button menu
function createMenu() {
    let tags = getUniqueTags();

    tags.forEach(item => {
        $(".filter-menu")
        .append("<div class='filter filter-btn' data-tag='" + item + "'>" + item + "</div>")
    });
}

// Populate the tab menu for tutorial groups
function populateTabs() {
    let groups = getUniqueGroups();

    groups.forEach(item => {
        $(".tab-menu")
        .append("<div class='tab' data-group='" + item + "'>" + item + "</div>")
    });
}

$(document).ready(function() {
    createMenu();
    populateTabs();
});

// Now add filtering functionality to the tag buttons

let selectedTagSet = new Set();

// Handle tag button selection and filtering
$(document).on("click", ".filter-btn", function () {
    console.log("button clicked!"); // for debugging
    let parent = $(this).closest(".filter-menu");
    let tag = $(this).data("tag");

    if (tag === "all") {
        if (!parent.hasClass("all-tag-selected")) {
            selectedTagSet.clear();
            parent.addClass("all-tag-selected");
            $(".filter-btn").removeClass("selected"); // Deselect all buttons
        }
    } else {
        if (selectedTagSet.has(tag)) {
            selectedTagSet.delete(tag);
            $(this).removeClass("selected"); // Deselect button

            if (selectedTagSet.size === 0) {
                parent.addClass("all-tag-selected");
            }
        } else {
            parent.removeClass("all-tag-selected");
            selectedTagSet.add(tag);
            $(this).addClass("selected"); // Highlight selected button
        }
    }

    filterCards(); // Trigger filtering immediately
});

// Function to filter cards based on both tags and groups
function filterCards() {
    $(".card").each(function () {
        let tagsData = $(this).data("tags") ?? "";
        let cardTags = tagsData.split(",").map(tag => tag.trim());
        let groupName = $(this).data("model-group-name");

        let matchesTags = selectedTagSet.size === 0 || [...selectedTagSet].every(tag => cardTags.includes(tag));
        let matchesGroup = selectedGroup === "all" || groupName === selectedGroup;

        console.log("Card tags:", cardTags);  // Log card tags
        console.log("Selected tags:", selectedTagSet);  // Log selected tags
        console.log("Matches tags:", matchesTags);  // Log matchesTags
        console.log("Matches group:", matchesGroup);  // Log matchesGroup

        $(this).toggleClass("hidden", !(matchesTags && matchesGroup));
    });
}

// Add similar filtering functionality to the model group tabs (but only single select)
let selectedGroup = "all"

// Handle tab selection and filtering
$(document).on("click", ".tab", function () {
    let group = $(this).data("group");

    if (group !== selectedGroup) {
        $(".tab").removeClass("tab-selected"); // deselect current tab
        $(this).addClass("tab-selected");
        selectedGroup = group
    }

    filterCards(); // Trigger filtering immediately
});
