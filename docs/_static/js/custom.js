// filtering process inspired by PyTorch Lightning sphinx theme

// Build an array from each tag that's present
var tagList = $(".tutorials-card-container").map(function() {
    return $(this).data("tags").split(",").map(function(item) {
        return item.trim();
    });
}).get();

function unique(value, index, self) {
    return self.indexOf(value) == index && value != "";
}

// Only return unique tags
var tags = tagList.sort().filter(unique);

// Add filter buttons to the top of the page for each tag
function createTagMenu() {
    tags.forEach(function(item){
        $(".tutorial-filter-menu").append(" <div class='tutorial-filter filter-btn filter' data-tag='" + item + "'>" + item + "</div>")
    });
}

createTagMenu();

$(document).on("click", ".page", function() {
    $('html, body').animate(
      {scrollTop: $("#dropdown-filter-tags").position().top},
      'slow'
    );
});
