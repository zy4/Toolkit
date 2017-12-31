/* Data table helpers */


/* REGEX FOR DB IDENTIFICATION*/
var uniprotReg = "^([A-Z0-9]{10}|[A-Z0-9]{6})$";
var scopReg = "^([defgh][0-9a-zA-Z\.\_]+)$";
var ecodReg = "^(ECOD_)";
var mmcifReg = "^(...._[a-zA-Z])$";
var mmcifShortReg = "^([0-9]+)$";
var pfamReg = "^(pfam[0-9]+&|^PF[0-9]+(\.[0-9]+)?)$";
var ncbiReg = "^([A-Z]{2}_?[0-9]+\.?\#?([0-9]+)?|[A-Z]{3}[0-9]{5}?\.[0-9])$";


function download(filename, text){
    var blob = new Blob([text], {type: "octet/stream"});
    if(window.navigator.msSaveOrOpenBlob) {
        window.navigator.msSaveBlob(blob, filename);
    } else {
        var a = document.createElement("a");
        a.style = "display: none";
        a.href = window.URL.createObjectURL(blob);
        a.download = filename;
        // Append anchor to body.
        document.body.appendChild(a);
        a.click();
        window.URL.revokeObjectURL(a.href);
        // Remove anchor from body
        a.remove();
    }
    $.LoadingOverlay("hide");
}



/* Slider */
function slider_show(sequence_length, start, end) {

    var tooltip = $('<div id="tooltip" />').css({
        position: 'absolute',
        top: -20
    }).show();

    var tooltip2 = $('<div id="tooltip2" />').css({
        position: 'absolute',
        top: -20

    }).show();



    $("#flat-slider").slider({
        range: true,
        orientation: 'horizontal',
        min: 1,
        max: sequence_length,
        step: 1,
        values: [start, end],
        slide: function(event, ui) {
            tooltip.text(ui.values[0]);
            tooltip2.text(ui.values[1]);
        },
        change: function(event, ui) {
            var sliderCoords =  $('#flat-slider').slider("option", "values");
        }
    });

    tooltip.text(start);
    tooltip2.text(end);

    $("#flat-slider").slider({}).find(".ui-slider-handle:first").append(tooltip);

    $("#flat-slider").slider({}).find(".ui-slider-handle:last").append(tooltip2);

}

function getSliderRange() {

    return $('#flat-slider').slider("option", "values");
}

function resubmitSection(sequence, name) {
    var sliderRange = getSliderRange();
    var resubmitSeqs = [];

    resubmitSeqs.push(name + '\n');
    resubmitSeqs.push(sequence.substr(sliderRange[0], sliderRange[1]) + '\n');

    $('a[href="#tabpanel-Input"]').click();
    $('#alignment').val(resubmitSeqs.join(''));
}


/**
 * forward the given forwardData to the given tool
 * @param tool
 * @param forwardData
 */
function forward(tool, forwardData){
    if(forwardData === ""){
        alert("No sequence(s) selected!");
        $.LoadingOverlay("hide");
        return;
    }
    try {
        localStorage.setItem("resultcookie", forwardData);
        window.location.href = "/#/tools/" + tool;
    } catch(e) {
        if (isQuotaExceeded(e)) {
            // Storage full, maybe notify user or do some clean-up
            $.LoadingOverlay("hide");
            alert("File is too big to be forwarded!" )
        }

    }
}

function forwardPath(tool, forwardPath){
    m.route("/tools/" + tool);
    $.ajax({
        type: 'GET',
        contentType: "charset=utf-8",
        url: forwardPath,
        error: function(){
            $.LoadingOverlay("hide");
        }
    }).done(function (data) {
    if(tool === "reformat"){
            setInterval(function(){ myCodeMirror.setValue(data); $.LoadingOverlay("hide"); }, 100);
        }
    else {
            $('#alignment').val(data);
        }
        validationProcess($('#alignment'),$("#toolnameAccess").val());
        $.LoadingOverlay("hide");
    })
}

function isQuotaExceeded(e) {
    var quotaExceeded = false;
    if (e) {
        if (e.code) {
            switch (e.code) {
                case 22:
                    quotaExceeded = true;
                    break;
                case 1014:
                    // Firefox
                    if (e.name === 'NS_ERROR_DOM_QUOTA_REACHED') {
                        quotaExceeded = true;
                    }
                    break;
            }
        } else if (e.number === -2147024882) {
            // Internet Explorer 8
            quotaExceeded = true;
        }
    }
    return quotaExceeded;
}
// load forwarded data into alignment field
$(document).ready(function() {
    var resultcookie = localStorage.getItem("resultcookie");
    $('#alignment').val(resultcookie);
    localStorage.removeItem("resultcookie");
    $.LoadingOverlay("hide");
});


function identifyDatabase(id){
    if (id === null)
        return null;
    if(id.match(scopReg))
        return "scop";
    else if(id.match(ecodReg))
        return "ecod";
    else if(id.match(mmcifShortReg))
        return "mmcif";
    else if(id.match(mmcifReg))
        return "mmcif";
    else if(id.match(pfamReg))
        return "pfam";
    else if(id.match(ncbiReg))
        return "ncbi";
    else if(id.match(uniprotReg))
        return "uniprot";
    else
        return null;
}




/* draw hits */
function drawHits(id, hits) {
    var s = Snap("#" + id);
    // Get maximum of query_end with map and reduce
    document.getElementById(id).setAttribute("viewBox", "0 0 " + hits.map(function (hit) {
            return hit.query_end;
        })
            .reduce(function (a, b) {
                return Math.max(a, b);
            }) + " " + hits.length);
    for (var i = 0; i < hits.length; ++i) {
        var hit = hits[i];
        var diff = hit.query_end - hit.query_begin;
        var r = s.rect(hit.query_begin, i, diff, 0.5, 0.5, 0.5);
        var text = s.text(diff / 2, i + 0.4, hit.struc);
        text.attr({
            'font-size': 0.5,
            'fill': 'white',
            'font-weight': 'bold'
        });
        var colors = calcColor(hit.prob).map(function (val) {
            return val * 255
        });
        r.attr({
            fill: Snap.rgb(colors[0], colors[1], colors[2])
        });
    }
}


function calcColor(prob) {
    var red, grn, blu;
    if (prob > 40) {
        var signif = ((prob - 40) / 60);
        var col = 4 * signif;
        if (col > 3) {
            col -= 3.0;
            red = 1;
            grn = 0.7 * (1 - col);
            blu = 0;
        } else if (col > 2) {
            col -= 2.0;
            red = col;
            grn = 0.7 + 0.3 * (1 - col);
            blu = 0;
        } else if (col > 1) {
            col -= 1.0;
            red = 0;
            grn = 1 - 0.3 * (1 - col);
            blu = 1 - col;
        } else {
            red = 0;
            grn = 0.7 * col;
            blu = 1;
        }
    } else {
        signif = Math.pow((prob / 40), 3);
        red = 0.2 * (1 - signif);
        grn = 0.2 * (1 - signif);
        blu = 0.2 + 0.8 * signif;
    }
    return [red, grn, blu];
}


function scrollToElem(num){
    num = parseInt(num);
    var elem = $('#tool-tabs').hasClass("fullscreen") ? '#tool-tabs' : 'html, body';
    if (num > shownHits) {
        $.LoadingOverlay("show");
        getHits(shownHits, num, wrapped, false).done(function(data){
            var pos = $('.aln"][value=' + num + ']').offset().top;
            $(elem).animate({
                scrollTop: pos - 100
            }, 1)
        }).then(function(){
            $.LoadingOverlay("hide");
        });
        shownHits = num;
    }else{
        var pos = $('.aln[value=' + num + ']').offset().top;
        $(elem).animate({
            scrollTop: pos - 100
        }, 1)
    }
}

function scrollToSection(name) {
    var elem = $('#tool-tabs').hasClass("fullscreen") ? '#tool-tabs' : 'html, body';
    var pos = $('#tool-tabs').hasClass("fullscreen") ? $('#' + name).offset().top + $(elem).scrollTop() : $('#' + name).offset().top + 25;
    $(elem).animate({
        scrollTop: pos
    }, 'fast');

}

function selectFromArray(checkboxes){
    _.range(1, numHits+1).forEach(function (currentVal) {
        $('input:checkbox[value='+currentVal+'][name="alignment_elem"]').prop('checked', checkboxes.indexOf(currentVal) !== -1 ? true : false);
    })
}

function hitlistBaseFunctions(){
    $(document).ready(function() {
        // add tooltipster to visualization
        $('#blastviz').find('area').tooltipster({
            theme: ['tooltipster-borderless', 'tooltipster-borderless-customized'],
            position: 'bottom',
            animation: 'fade',
            contentAsHTML: true,
            debug: false,
            maxWidth: $(this).innerWidth() * 0.6
        });
        $.LoadingOverlay("hide");
        followScroll(document);

        // add slider val
        $('.slider').on('moved.zf.slider', function () {
            $('#lefthandle').html($('#hidden1').val()).css({
                'color': 'white',
                'font-weight': 'bold',
                'padding-left': '2px'
            });
            $('#righthandle').html($('#hidden2').val()).css({
                'color': 'white',
                'font-weight': 'bold',
                'padding-left': '2px'
            });
        });
    });
}


Array.prototype.removeDuplicates = function () {
    return this.filter(function (item, index, self) {
        return self.indexOf(item) === index;
    });
};


function getsHitsManually(){
    if (!loading) {
        var end = shownHits + 100;
        end = end < numHits ? end : numHits;
        if (shownHits !== end) {
            getHits(shownHits, end, wrapped, false);
        }
        shownHits = end;
    }
}

function linkCheckboxes(){
    $('input:checkbox').on('change',function (e) {
        var currentVal = $(this).val();
        var currentState = $(this).prop('checked');

        // link checkboxes with same value
        $('input:checkbox[value=' + currentVal + '][name=alignment_elem]').each(function () {
            $(this).prop('checked', currentState);
        });

        if (currentState) {
            // push num of checked checkbox into array
            checkboxes.push(parseInt(currentVal));
            // make sure array contains no duplicates
            checkboxes = checkboxes.filter(function (value, index, array) {
                return array.indexOf(value) === index;
            });
        } else {
            // delete num of unchecked checkbox from array
            checkboxes = checkboxes.filter(function(x){return x !== currentVal});
        }

    });
}


function generateFilename(){
    return Math.floor(100000 + Math.random() * 900000).toString();
}

/**
 * wraps sequences for search tools
 * for this it empties the table "#alignmentTable"
 * and calls get Hits taking the boolean wrapped as a parameter
 */
function wrap(){
    wrapped = !wrapped;
    var elemArr =  $(".aln").toArray();
    var num = 1;
    for(var i =0 ; i < elemArr.length; i++){
        if($(elemArr[i]).isOnScreen()){
            num  = $(elemArr[i]).attr("value");
            break;
        }
    }
    $("#wrap").toggleClass("colorToggleBar");
    $("#wrap").toggleText("Unwrap Seqs", "Wrap Seqs");
    $("#alignmentTable").empty();
    getHits(0, shownHits, wrapped, false).then(function(){
        linkCheckboxes();
        scrollToElem(num);
    });

}

function colorAA(){
    colorAAs = !colorAAs;
    $.LoadingOverlay("show");
    $(".colorAA").toggleClass("colorToggleBar");
    var elemArr =  $(".aln").toArray();
    var num = 1;
    for(var i =0 ; i < elemArr.length; i++){
        if($(elemArr[i]).isOnScreen()){
            num  = $(elemArr[i]).attr("value");
            break;
        }
    }
    $("#alignmentTable").empty();
    getHits(0, shownHits, wrapped, colorAAs).then(function(){
        $.LoadingOverlay("hide");
        linkCheckboxes();
        scrollToElem(num);
    });
}

$.fn.isOnScreen = function(){
    var viewport = {};
    viewport.top = $(window).scrollTop();
    viewport.bottom = viewport.top + $(window).height();
    var bounds = {};
    bounds.top = this.offset().top;
    bounds.bottom = bounds.top + this.outerHeight();
    return ((bounds.top <= viewport.bottom) && (bounds.bottom >= viewport.top));
};
