// to do
//   -boxplot
//   -starplot
//   -subDir
//   -colnames
//   -check sorting

/* Table initialisation */
$(document).ready(function() {
	$('#main').dataTable( {
		"bProcessing": true,
        "aaData": data,
        "aoColumns": [
            { "mData": "gId",
              "mRender": function ( data, type, row ) {
                  href="http://www.ensembl.org/Homo_sapiens/Gene/Summary?g=" + data;
                  return '<a href="' + href + '">' + data + '</a>';
            }},            
            { "mData": "gName" },
            { "mData": "nOfT" },
            { "mData": "C1_tId" },
            { "mData": "C1_principal",
              "mRender": function ( data, type, row ) {
                  href="http://appris.bioinfo.cnio.es/report.html?id=" + row.C1_tId + 
                    "&namespace=Ensembl_Gene_Id&specie=homo_sapiens";
                  return '<a href="' + href + '">' + data + '</a>';
            }},
            { "mData": "C1_biotype" },
            { "mData": "C1_tExp" },
            { "mData": "C1_gExp" },
            { "mData": "C1_breadth" },

			{ "mData": "C2_tId" },
            { "mData": "C2_principal",
              "mRender": function ( data, type, row ) {
                  href="http://appris.bioinfo.cnio.es/report.html?id=" + row.C2_tId + 
                    "&namespace=Ensembl_Gene_Id&specie=homo_sapiens";
                  return '<a href="' + href + '">' + data + '</a>';
            }},
            { "mData": "C2_biotype" },
            { "mData": "C2_tExp" },
            { "mData": "C2_gExp" },
            { "mData": "C2_breadth" },
    
            { "mData": "pIdentity",
              "mRender": function ( data, type, row ) {
                  href="fixme" + row.gId+ 
                    "fixme";
                  // my $url="$info{'out_dir'}/prot_aln/$subDir/$gId.needle_mod.out";
                  return '<a href="' + href + '">' + data + '</a>';
            }},
            { "mData": "pdbEntry",
              "mRender": function ( data, type, row ) {
                  if (data != "NO") {
                        href="http://www.ebi.ac.uk/pdbe/widgets/unipdb?uniprot=" + data;
                        return '<a href="' + href + '">YES</a>';
                  } else {
                        return data;
                  }
            }},
            // { "mData": "boxplot" },
            // { "mData": "starplot" },
            { "mData": "rank" }
        ],
		"sDom": "<'row-fluid'<'span6'T><'span6'f>r>t<'row-fluid'<'span6'i><'span6'p>>",
		"oTableTools": {
			"aButtons": [ "print" ]
		}
	} );
} );

// $(document).ready(function() {
// 	//data
//   $('#main').dataTable({
//   		"bProcessing": true,
// //         "aaData": data,
//         "aoColumns": [
//             { "mData": "render_engine" },
//             { "mData": "browser" },
//             { "mData": "platform" },
//             { "mData": "enging_version" },
//             { "mData": "css_grade" }
//         ]
//   });
// } );