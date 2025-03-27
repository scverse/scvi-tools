console.log("check card containers:", $(".card-container").length);
console.log("check filter menus:", $(".filter-menu").length);

// Create a set of all unique tags that appear in the cards
function getUniqueTags() {
    let tagSet = new Set();
    console.log("getUniqueTags is running"); // Debugging log

    $(".card-container").each(function() {
        let tags = $(this).data('tags').split(",")
        .map(tag => tag.trim()) // Get rid of whitespace
        .filter(tag => tag !== ""); // Make sure there are no empty tags

        console.log("Tags found:", tags); // See if tags are being detected

        tags.forEach(tag => tagSet.add(tag)); // Add tag to the set
    });

    let uniqueTags = Array.from(tagSet);
    uniqueTags.sort();
    return uniqueTags;
}

// Create the tag button menu
function createMenu() {
    let tags = getUniqueTags();

    tags.forEach(item => {
        $(".filter-menu")
        .append("<div class='filter filter-btn' data-tag='" + item + "'>" + item + "</div>")
    });
}

$(document).ready(function() {
    createMenu();
});

// Now add filtering functionality to the tag buttons

let selectedTagSet = new Set();

// Handle tag button selection and filtering
$(".filter-btn").on("click", function () {
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

// Function to filter cards
function filterCards() {
    if ($(".filter-menu").hasClass("all-tag-selected")) {
        $(".card").removeClass("hidden"); // Show all cards
    } else {
        $(".card").each(function () {
            let cardTags = $(this).data("tags").split(",").map(tag => tag.trim());

            // Show only if the card contains all selected tags
            let shouldShow = [...selectedTagSet].every(tag => cardTags.includes(tag));

            $(this).toggleClass("hidden", !shouldShow);
        });
    }
}
