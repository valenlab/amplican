// When click on the hide button on a header, hide the body
function hide(target) {

	console.log("Hide")
	
	if (target.style.display == 'none'){
		
		target.style.display = 'flex';

	}
            
    else{
		
        target.style.display = 'none';
        
	}

}
