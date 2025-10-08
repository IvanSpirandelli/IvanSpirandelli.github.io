document.addEventListener("DOMContentLoaded", function() {
    var coll = document.getElementsByClassName("collapsible");
    var i;

    for (i = 0; i < coll.length; i++) {
        coll[i].addEventListener("click", function() {
            this.classList.toggle("active");
            
            // --- CHANGED PART ---
            var content = this.nextElementSibling;
            // Loop until we find the 'content' div
            while (content && !content.classList.contains('content')) {
                content = content.nextElementSibling;
            }
            // --- END OF CHANGE ---

            if (content) { // Make sure we found the content div
                if (content.style.maxHeight){
                    content.style.maxHeight = null;
                } else {
                    content.style.maxHeight = content.scrollHeight + "px";
                } 
            }
        });
    }
});