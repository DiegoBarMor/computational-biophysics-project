document.addEventListener("DOMContentLoaded", function() {
    var selectElement = document.getElementById("page-select");
    var currentPage = window.location.pathname.split("/").pop().replace(".html", "");

    // Set the dropdown to the current page
    for (var i = 0; i < selectElement.options.length; i++) {
        if (selectElement.options[i].value == currentPage) {
            selectElement.selectedIndex = i;
            selectElement.title = selectElement.options[i].text;
            break;
        }
    }

    selectElement.addEventListener("change", function() {
        var selectedValue = selectElement.value;
        if (selectedValue == currentPage) return;
        selectElement.title = selectElement.options[selectElement.selectedIndex].text;
        window.location.href = "https://diegobarmor.github.io/computational-biophysics-project/" + selectedValue + ".html";
    });
});
