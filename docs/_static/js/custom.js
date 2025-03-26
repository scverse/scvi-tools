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
