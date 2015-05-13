
// http://127.0.0.1:7134/?adjustmethod=meancentring&batcheffect1=0&batcheffect2=0&batcheffect3=0&countA1=10&countA2=0&countA3=0&countB1=1&countB2=1&countB3=0&countC1=0&countC2=10&countC3=0&groupeffectA=0&groupeffectB=0&groupeffectC=0&indexgene=1&pair=AB&rngseed=139&

// from http://stackoverflow.com/questions/901115/how-can-i-get-query-string-values-in-javascript
var urlParams;
(window.onpopstate = function () {
    var match,
        pl     = /\+/g,  // Regex for replacing addition symbol with a space
        search = /([^&=]+)=?([^&]*)/g,
        decode = function (s) { return decodeURIComponent(s.replace(pl, " ")); },
        query  = window.location.search.substring(1);
				//alert("query="+query);
    urlParams = {};
    while (match = search.exec(query))
       urlParams[decode(match[1])] = decode(match[2]);
    
})();


function initparams()
{
  //alert("initparams");
  changevalue(urlParams);
  
  // Some settings only availible when app is run locally
  var host = window.location.hostname;
  if(host=="127.0.0.1" || host=="localhost")
  {
    document.getElementById('localhostonly').style.display = 'block'
  }
  setbatchrows();  
}

function setbatchrows()
{
  //only show a certain number of batch-controls
  var elm=document.getElementById('numberofbatches');
  visiblebatches=elm.options[elm.selectedIndex].value;
  for(var i=1;i<21;i++)
  {    
    if(i<=visiblebatches)
    {
      document.getElementById('batchrow'+i).style.display = 'block';
      document.getElementById('batcheffectvalue'+i).style.display = 'block';
      
    }
    else
    {
      document.getElementById('batchrow'+i).style.display = 'none';
      document.getElementById('batcheffectvalue'+i).style.display = 'none';
    }
  }
  
}
 
 function test1()
 {
 		alert("urlParams="+JSON.stringify(urlParams))
 }
 
 
 
 function changevalue(ar)
 {
  // console.log("i changevalue");
	for (var key in ar)
 	{     
     var elm = document.getElementById(key);
     //console.log("changem key="+key+"  val=" + ar[key]+ " typee="+ elm.type);
     if(elm)
     {
     		//alert("elmtype="+elm.type);
       if(elm.type=="checkbox")
       {
         //console.log("setter checkbox");
          elm.checked=false;
          if(ar[key]=="TRUE")
            elm.checked=true;
       }
       else
       {       
        document.getElementById(key).value=ar[key];	 	    
       }
       $('#'+key).trigger('change');
     }
  }


 }