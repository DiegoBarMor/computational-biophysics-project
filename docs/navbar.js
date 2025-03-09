document.addEventListener("DOMContentLoaded", function() {
    var selectElement = document.getElementById("page-select");
    selectElement.addEventListener("change", function() {
        window.location.href = "https://diegobarmor.github.io/computational-biophysics-project/" + selectElement.value + ".html";
    });
});
