var filterCategory = "";
var myTable = {
  "sDom": "<'row-fluid'<'span6'T><'span6'f>r>t<'row-fluid'<'span6'i><'span6'p>>",
  "oTableTools": {
    "aButtons": [ "print" ]
  },
  "bProcessing": true,
  "aaData": data,
  "aoColumns": [

      // general info
      { "sTitle": "gId",
        "mData": "gId",
        "mRender": function ( data, type, row ) {
            href="http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=" + data;
            return '<a href="' + href + '">' + data + '</a>';
      }},            
      { "sTitle": "gName",
        "mData": "gName" },
      { "sTitle": "nOfT",
        "mData": "nOfT" },

      // C1 info
      { "sTitle": "C1.tId",
        "mData": "C1_tId" },
      { "sTitle": "C1.principal",
        "mData": "C1_principal",
        "mRender": function ( data, type, row ) {
            href="http://appris.bioinfo.cnio.es/report.html?id=" + row.C1_tId + 
              "&namespace=Ensembl_Gene_Id&specie=homo_sapiens";
            return '<a href="' + href + '">' + data + '</a>';
      }},
      { "sTitle": "C1.biotype",
        "mData": "C1_biotype" },
      { "sTitle": "C1.tExp",
        "mData": "C1_tExp" },
      { "sTitle": "C1.gExp",
        "mData": "C1_gExp" },
      { "sTitle": "C1.breadth",
        "mData": "C1_breadth" },

      // C2 info
      { "sTitle": "C2.tId",
        "mData": "C2_tId" },
      { "sTitle": "C2.principal",
        "mData": "C2_principal",
        "mRender": function ( data, type, row ) {
            href="http://appris.bioinfo.cnio.es/report.html?id=" + row.C2_tId + 
              "&namespace=Ensembl_Gene_Id&specie=homo_sapiens";
            return '<a href="' + href + '">' + data + '</a>';
      }},
      { "sTitle": "C2.biotype",
        "mData": "C2_biotype" },
      { "sTitle": "C2.tExp",
        "mData": "C2_tExp" },
      { "sTitle": "C2.gExp",
        "mData": "C2_gExp" },
      { "sTitle": "C2.breadth",
        "mData": "C2_breadth" },

      // general info again
      { "sTitle": "pIdentity",
        "mData": "pIdentity",
        "mRender": function ( data, type, row ) {
          if (data != "NA") {
            string=row.gId;
            subdir=string.substring(0, 12);
            href="./data/prot_aln/" + subdir + "/" + row.gId + ".needle_mod.out";
            return '<a href="' + href + '">' + data + '</a>';
          } else {
            return data;
          }
      }},
      { "sTitle": "pdbEntry",
        "mData": "pdbEntry",
        "mRender": function ( data, type, row ) {
            if (data != "NO") {
                  href="http://www.ebi.ac.uk/pdbe/widgets/unipdb?uniprot=" + data;
                  return '<a href="' + href + '">YES</a>';
            } else {
                  return data;
            }
      }},
      { "sTitle": "distrplot",
        "mData": null,
        "mRender": function ( data, type, row ) {
            string=row.gId;
            subdir=string.substring(0, 12);
            href="./data/plots/distrplots/" + subdir + "/" + row.gId + ".pdf";
            return '<a href="' + href + '">distrplot</a>';
      }},
      { "sTitle": "starplot",
        "mData": null,
        "mRender": function ( data, type, row ) {
            string=row.gId;
            subdir=string.substring(0, 12);
            href="./data/plots/starplots/" + subdir + "/" + row.gId + ".pdf";
            return '<a href="' + href + '">starplot</a>';
      }},            
      { "sTitle": "rank",
        "mData": "rank" }
  ]
};

function filterTable (C1_biotype, C1_filter, C2_biotype, C2_filter) {
      if ( C1_biotype == C1_filter && C2_biotype == C2_filter ) {
        return true;
      }
      return false;
};

function filterTableOther (C1_biotype, C2_biotype) {
      if ( !(C1_biotype == "protein_coding" && C2_biotype == "protein_coding")
        && !(C1_biotype == "protein_coding" && C2_biotype == "nonsense_mediated_decay")
        && !(C1_biotype == "protein_coding" && C2_biotype == "retained_intron") 
        && !(C1_biotype == "protein_coding" && C2_biotype == "processed_transcript")
        && !(C1_biotype == "nonsense_mediated_decay" && C2_biotype == "protein_coding")
        && !(C1_biotype == "retained_intron" && C2_biotype == "protein_coding")
        && !(C1_biotype == "processed_transcript" && C2_biotype == "protein_coding") ) {
        return true;
      }
      return false;
};

$.fn.dataTableExt.afnFiltering.push(    
    function ( oSettings, aData, iDataIndex ) {
      // return all rows
      if (filterCategory === "") { return true };

      // filtering
      var C1_biotype = aData[5];
      var C2_biotype = aData[11];
      var C1_filter;
      var C2_filter;

      if (filterCategory == "pc_to_pc") {
        C1_filter="protein_coding";
        C2_filter="protein_coding";
        return filterTable(C1_biotype, C1_filter, C2_biotype, C2_filter);
      }
      if (filterCategory == "pc_to_nmd") {
        C1_filter="protein_coding";
        C2_filter="nonsense_mediated_decay";
        return filterTable(C1_biotype, C1_filter, C2_biotype, C2_filter);
      }
      if (filterCategory == "pc_to_ri") {
        C1_filter="protein_coding";
        C2_filter="retained_intron";
        return filterTable(C1_biotype, C1_filter, C2_biotype, C2_filter);
      }
      if (filterCategory == "pc_to_pt") {
        C1_filter="protein_coding";
        C2_filter="processed_transcript";
        return filterTable(C1_biotype, C1_filter, C2_biotype, C2_filter);
      }   
      if (filterCategory == "nmd_to_pc") {
        C1_filter="nonsense_mediated_decay";
        C2_filter="protein_coding";
        return filterTable(C1_biotype, C1_filter, C2_biotype, C2_filter);
      }               
      if (filterCategory == "ri_to_pc") {
        C1_filter="retained_intron";
        C2_filter="protein_coding";
        return filterTable(C1_biotype, C1_filter, C2_biotype, C2_filter);
      }   
      if (filterCategory == "pt_to_pc") {
        C1_filter="processed_transcript";
        C2_filter="protein_coding";
        return filterTable(C1_biotype, C1_filter, C2_biotype, C2_filter);
      }   
      if (filterCategory == "other") {
        return filterTableOther(C1_biotype, C2_biotype);
      }
    }
  );

function updateTable ( filter ) {
  filterCategory = filter;
  var oTable = $('#main').dataTable();
  oTable.fnDraw();
};

// Table initialisation
$(document).ready(function() {
	$('#main').dataTable( myTable );
} );