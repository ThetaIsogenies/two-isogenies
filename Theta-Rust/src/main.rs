#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

fn sqi_example() {
    use theta_rs::theta254::{product_isogeny, CouplePoint, Curve, EllipticProduct, Fq, Point};

    println!("Testing chain for 254-bit parameters");

    // Curve coefficients
    let A1 = Fq::ZERO;
    let A2_str = "2687f041b47d1ea9c00ae938b2a761aac6bad9ea80dd9dbe1c24d9ef697d7d0475898661998dd3a7b186e2558d1cf0dd771fb49483988c2ff578547815e8f00e";
    let (A2, _) = Fq::decode(&hex::decode(A2_str).unwrap());

    // Kernel Data
    let P1_X_str = "63385515142ea2f7a0c3747c22668aef99b31f43987354116da1915c24ff2f0a279051962feace976834986fc955b11bb8e1ffea47d8ce994ad22ee86c7f7a00";
    let P1_Y_str = "38243631804d307b72ecd037da591d1f06ca606bf1bb71d77ce10467b00f7b082472068bfea9eb9d80b68e04bb194a23e5214ba41625915d8e590024e5dcf611";
    let P2_X_str = "3e218d8b18cf09ce29cefaa35467225134910411fe33625136e50f7b9c59e51c517629d8786a98603cc06470a10dea83f7eca03b9b378297c21755bf0aee1324";
    let P2_Y_str = "06d56830cea82c91cf4059078566145d4b90d992177916ffe1380060057d75277610f27fc5558e5d028699493a300d84521f0c077c6e52c6adc1820c8f53f11a";
    let Q1_X_str = "3b569c320301eb5aa1aa078b7399d31e2e0e6c70e91223d1d3346be7145c100b16e4c042048452157133174122c04b1c9a17f38c28e959828933a95eebf6a305";
    let Q1_Y_str = "04a2a994555e67a3ddbb5f87dbef7903a9fc9724cb36e51924d28522222c3d2bd07a4d44b06b176b63d78733a1a4839606276b8c523bd6dc8f23de01e2817e17";
    let Q2_X_str = "4076f17bf841c5d6acc194e75a4cd020fe3ea03b0914cf3f3db9cc882f8b4724a18ce4bb13c2b3c46abd8e6cdb502dd7f48e58ee49d6c1d632532f6ed995e12b";
    let Q2_Y_str = "5b88c588b1f27496800acbef34817c5fa5cbcd728de00e31a46fc7aa6ef9af0e173d1ba96e1b2c2ebc5bc3dd3f980344b508b9df1863fb624855dc1a8cc17b12";
    let (P1_X, _) = Fq::decode(&hex::decode(P1_X_str).unwrap());
    let (P1_Y, _) = Fq::decode(&hex::decode(P1_Y_str).unwrap());
    let (P2_X, _) = Fq::decode(&hex::decode(P2_X_str).unwrap());
    let (P2_Y, _) = Fq::decode(&hex::decode(P2_Y_str).unwrap());
    let (Q1_X, _) = Fq::decode(&hex::decode(Q1_X_str).unwrap());
    let (Q1_Y, _) = Fq::decode(&hex::decode(Q1_Y_str).unwrap());
    let (Q2_X, _) = Fq::decode(&hex::decode(Q2_X_str).unwrap());
    let (Q2_Y, _) = Fq::decode(&hex::decode(Q2_Y_str).unwrap());

    // Points to push through isogeny
    let PA_X_str = "95a2c9abfff13d34db8c9d79ecdeb8cf8cbdb2de119f5fc3c51e69dcffab602e79fc55b299bdc279f8c3c20eac06322b43cb71718e367e1f795b4c8fc59fbb0b";
    let PA_Y_str = "b5b5fdfca0cc051994cbeb4338a5e629714df6cc7496cdd98900fdf281ff342990f21fa62afcf762a1c6f3e65635ab87f5c90269722434f2479482b2ec089a05";
    let (PA_X, _) = Fq::decode(&hex::decode(PA_X_str).unwrap());
    let (PA_Y, _) = Fq::decode(&hex::decode(PA_Y_str).unwrap());

    // Result to check against (Debugging)
    let im_PA_X_str = "0b2588ba81d363ee3d7c11e52e8e9db94986820a433eefb9c294de6393a1fa09eaa2baf6945de410af1ab4c189a1d5b885257ad7d2133400b03a24e4aef1f603";
    let im_PA_Y_str = "bc5f46a11f3497cdbf9c05497ad8fa57104bdfb9431c3c327b702d164618d8201cf44b3eebc8a5d7c7214668523ec2e8fd4ad135faed66074e605abb098d1409";
    let im_PB_X_str = "28f000f5b02a6d3e96bc1e53c3d3eb63c77d9ac8b01386da6a1ff31dff0df101103ce4b37cccfd5a9d6a67fb207e69d12317c5eb2a94a09ae156f343c6c5f30d";
    let im_PB_Y_str = "12a4e94423744ae02c830a5188e345b9eac47823cfd7a70062995a0517137a028863d606013f8d4c42018cd7b540ad5b8a7135bbbfa3caf94019832346732010";
    let (im_PA_X, _) = Fq::decode(&hex::decode(im_PA_X_str).unwrap());
    let (im_PA_Y, _) = Fq::decode(&hex::decode(im_PA_Y_str).unwrap());
    let (im_PB_X, _) = Fq::decode(&hex::decode(im_PB_X_str).unwrap());
    let (im_PB_Y, _) = Fq::decode(&hex::decode(im_PB_Y_str).unwrap());

    // Curves which define elliptic product
    let E1 = Curve::new(&A1);
    let E2 = Curve::new(&A2);
    let E1E2 = EllipticProduct::new(&E1, &E2);

    // Points on E1 x E2
    let P1 = Point::new_xy(&P1_X, &P1_Y);
    let P2 = Point::new_xy(&P2_X, &P2_Y);
    let Q1 = Point::new_xy(&Q1_X, &Q1_Y);
    let Q2 = Point::new_xy(&Q2_X, &Q2_Y);
    let P1P2 = CouplePoint::new(&P1, &P2);
    let Q1Q2 = CouplePoint::new(&Q1, &Q2);

    // Point to push through isogeny
    let PA = Point::INFINITY;
    let PB = Point::new_xy(&PA_X, &PA_Y);
    let PAPB = CouplePoint::new(&PA, &PB);
    let image_points = [PAPB];

    // Points to compare against
    let im_PA = Point::new_xy(&im_PA_X, &im_PA_Y);
    let im_PB = Point::new_xy(&im_PB_X, &im_PB_Y);

    // Length of isogeny chain
    let n = 126;

    // Strategy computed from strategy.py
    let strategy: [usize; 125] = [
        55, 34, 21, 15, 5, 3, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1,
        1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1,
        5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1,
        1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3,
        2, 1, 1, 1, 1, 1,
    ];

    // Compute chain
    let (E3E4, images) = product_isogeny(&E1E2, &P1P2, &Q1Q2, &image_points, n, &strategy);

    let (E3, E4) = E3E4.curves();
    println!("E3: {}", E3);
    println!("E4: {}", E4);

    let (P1, P2) = images[0].points();

    let P1_check: bool = P1.equals(&im_PA) | P1.equals(&(-im_PA)) == u32::MAX;
    let P2_check: bool = P2.equals(&im_PB) | P2.equals(&(-im_PB)) == u32::MAX;

    println!("First image matches what is expected: {}", P1_check);
    println!("Second image matches what is expected: {}", P2_check);
    println!();
}

fn three_lambda_example() {
    println!("Testing chain for 381-bit parameters");

    use theta_rs::theta381::{product_isogeny, CouplePoint, Curve, EllipticProduct, Fq, Point};
    // Curve coefficients
    let A1 = Fq::ZERO;
    let A2_str = "346bc58bbc832c1c6c4c8b6cfeb58a483385fcb98594ac4d22e85ce21a520f12138a05fa5cb716068ef218ee1fefb012476841a87b5755af5d1c7041ce205487d5a98922d6181f68b8ec6a41058b8501baab5bc935017509142d7e1a212b2209";
    let (A2, _) = Fq::decode(&hex::decode(A2_str).unwrap());

    // Kernel Data
    let P1_X_str = "3a1111b79a2ac3a148c08c10cfe6f67ad7a9bd16aa56458013a99f266cbe4146e183bbc865881aad06022551b0e87f0aac35b5ca7f605e29cf7fa11f130a77a63ccc7a3930ffb0624af8856c829f90f70b21a87295b14b7e3b8cec6b23fd3f06";
    let P1_Y_str = "8cd2e58f9fa3f3678914731daf18d4fc33ac1fc491dfcd2bb0f51f4029627bd4e324205cda127a52a4fbab7d3979920b215ad30c48034e000579fecc8e26ec908c4348c98b0a2fbce0a54556e6d3489db11bfd7e2804077e18448f7dc945d008";
    let P2_X_str = "c055d51a89fbd10363d07613b1a71c925750811a11207d583642b71dd296c217e576205442b87c7d0ecf00f29d67e90b8f6898353649e1d53cd75601a1a5148384fa2a9fc89084746fc0d10b68e3d74e0f1ec66c74a582f98ffa8ba2497e080b";
    let P2_Y_str = "542497699a437f71c004bb51a0033e5199257811be7b397b84a6195532f8ec90de59953a1cd9b17bad40fa5c7897e5086ca8159f41a867455f53a1ea6766319bd3286eae4c30c37c1bdedbaa873ab854962f6526f82469b23d73194e1d9f9908";
    let Q1_X_str = "aed83312515b93af8e0fda085d3f7f96c455fe3a7a60e8c45239abb008f9c483eb61225a4ff0410be385c559766df60e28bae5c902ee1c46685717687ab6b5d67e5a87402092578d47842cd2ff943104ca7adff59544f7124d011919c564d607";
    let Q1_Y_str = "f8edffdaa53ac760c85b496ed64acdda15aed38d4b1e729f171caa77f482c406db35ac53cb94b8f9631d66c7b1b75310735f6a1d09815dadf8a1f849cdbee90f8d0cd0d2ef77f1a4b7ea303f00454f7b718e90433a487bc1f7e75d1d63106911";
    let Q2_X_str = "d9e22eb75836ea5dd0e5f08cc939b63976c948f6c1bf71b3568598c7c752df2a4ea0768d978d62895c584b1588fc9c111f36bacc042a17fea60fdfeb10608120a7afad2f61064b7fb7dd4e688089348f00d51e5182c51a88bfeff67f68564909";
    let Q2_Y_str = "79c00be165d30dcc11c5a46f512cdec224a09dc3213a0834318aba1ca33e99d4e1b992e91ec91faf714a85e7a1cfe9138e36f63a53f494223dc07a38f79c9edb090806959f0ad779b432768a9274ede5db92390a9040fd1b57e78baee33fc405";
    let (P1_X, _) = Fq::decode(&hex::decode(P1_X_str).unwrap());
    let (P1_Y, _) = Fq::decode(&hex::decode(P1_Y_str).unwrap());
    let (P2_X, _) = Fq::decode(&hex::decode(P2_X_str).unwrap());
    let (P2_Y, _) = Fq::decode(&hex::decode(P2_Y_str).unwrap());
    let (Q1_X, _) = Fq::decode(&hex::decode(Q1_X_str).unwrap());
    let (Q1_Y, _) = Fq::decode(&hex::decode(Q1_Y_str).unwrap());
    let (Q2_X, _) = Fq::decode(&hex::decode(Q2_X_str).unwrap());
    let (Q2_Y, _) = Fq::decode(&hex::decode(Q2_Y_str).unwrap());

    // Points to push through isogeny
    let PA_X_str = "49c59ebe43e505e387d8306fcba301dbec0499bb5a0518e06c1feb4cabf78a8d5588089f6f7f724776e484efd54f5f0270399cd2318d17876bd7a6e0b75891d0b2b4489a85af988b7518f8b05bb8bdeb2fbc2433dfb739861ca2471fccef4d0c";
    let PA_Y_str = "47d0dde303fbc0e86869647e6c28e82bb82250761f5dce33bc4bb1a895721832270e7a063d98fb1eb1bade0760cf3811269db26ba866299486ae92a8182afcccbacaa6bce40d2905949ff30ef598e1f2a3f1ba01aad7aaac98f256263da0c70b";
    let (PA_X, _) = Fq::decode(&hex::decode(PA_X_str).unwrap());
    let (PA_Y, _) = Fq::decode(&hex::decode(PA_Y_str).unwrap());

    // Result to check against (Debugging)
    let im_PA_X_str = "530fc1ec207f15cb6fa60b1b2d242806fb743898d1267644c474c6321d41412d163209f8b9402fa8b567e6638bce1e0679bfbd40d3451406b5ef47652a162d567bcce1fbec10db716ef3db52dbdd466c78d3fcf93475307072a59a894acb1811";
    let im_PA_Y_str = "32f890c6a1c5d88ed1eb0c51aa3a9c68b29b97b6747161db58d264dffef926428ebdebe0bdd0b612b8ad47f23f189a0950bc7f185216d104f8fe601eefd68089edfc217ee980e60a6286d92b23fa08f47e72383cc19fbf91e1fdfa49065b4504";
    let im_PB_X_str = "a7de3a5eeb67be3f855942987f488c5ab813a083332846897741870a76356a12cc5f765774ebba539681369b90cabc01f6e59bcb0a300650b9b0865484f9ad735a55b7fa098c970febef388e05f1af9c2217278f6e3c6124178b79215e78c705";
    let im_PB_Y_str = "f2e6cb08c0a84323c8983727cacebfb3b73bb999507eb1a5f55a4519f02dee4f2362c2da220bb81ac010a59df128471092696fe93e5b9bceaede6dc2825f7d8920dc807b5529f4a65aa4f6179a8dc60adc3b461dc464f73c65d41b71267e150b";
    let (im_PA_X, _) = Fq::decode(&hex::decode(im_PA_X_str).unwrap());
    let (im_PA_Y, _) = Fq::decode(&hex::decode(im_PA_Y_str).unwrap());
    let (im_PB_X, _) = Fq::decode(&hex::decode(im_PB_X_str).unwrap());
    let (im_PB_Y, _) = Fq::decode(&hex::decode(im_PB_Y_str).unwrap());

    // Curves which define elliptic product
    let E1 = Curve::new(&A1);
    let E2 = Curve::new(&A2);
    let E1E2 = EllipticProduct::new(&E1, &E2);

    // Points on E1 x E2
    let P1 = Point::new_xy(&P1_X, &P1_Y);
    let P2 = Point::new_xy(&P2_X, &P2_Y);
    let Q1 = Point::new_xy(&Q1_X, &Q1_Y);
    let Q2 = Point::new_xy(&Q2_X, &Q2_Y);
    let P1P2 = CouplePoint::new(&P1, &P2);
    let Q1Q2 = CouplePoint::new(&Q1, &Q2);

    // Point to push through isogeny
    let PA = Point::INFINITY;
    let PB = Point::new_xy(&PA_X, &PA_Y);
    let PAPB = CouplePoint::new(&PA, &PB);
    let image_points = [PAPB];

    // Points to compare against
    let im_PA = Point::new_xy(&im_PA_X, &im_PA_Y);
    let im_PB = Point::new_xy(&im_PB_X, &im_PB_Y);

    // Length of isogeny chain
    let n = 208;

    // Strategy computed from strategy.py
    let strategy: [usize; 207] = [
        84, 55, 34, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5,
        3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1,
        1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2,
        1, 1, 1, 1, 1, 29, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5,
        3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1,
        1, 8, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1,
    ];

    // Compute chain
    // let (E3E4, images) = product_isogeny(&E1E2, &P1P2, &Q1Q2, &mut image_points, n, &strategy);
    let (E3E4, images) = product_isogeny(&E1E2, &P1P2, &Q1Q2, &image_points, n, &strategy);

    let (E3, E4) = E3E4.curves();
    println!("E3: {}", E3);
    println!("E4: {}", E4);

    let (P1, P2) = images[0].points();

    let P1_check: bool = P1.equals(&im_PA) | P1.equals(&(-im_PA)) == u32::MAX;
    let P2_check: bool = P2.equals(&im_PB) | P2.equals(&(-im_PB)) == u32::MAX;

    println!("First image matches what is expected: {}", P1_check);
    println!("Second image matches what is expected: {}", P2_check);
    println!();
}

fn festa_example() {
    println!("Testing chain for FESTA parameters");

    use theta_rs::thetaFESTA::{product_isogeny, CouplePoint, Curve, EllipticProduct, Fq, Point};
    // Curve coefficients
    let A1_str = "d268d78b4a92566637452bf184ade847911397b523d84750d9beddf94d402fa5ca047b8a729844f8362094d212cd6fb9121a44241fbc8a0cc2bac1562744f015e3d41a6fa75014fa6a00ad4ac9ea34b202160d0876574324a65b7b2253142a6f408e89094e8ed1ba44fb56eae37cd6de728f177c83391c6763bbf493fc33874cee9259989d59d403886d48dad157f40238acb776ecf4e38fb5653ddada4afd6d0402099fc10c379b5190c8c195531815b4729da27d43cbaf41b5d53185cabaff515669eb38ef7a39eb1ecbae997f063f47fa9952aef62eed5ae8b647d578dcc1da25aaf197d1c2cfcb9043e9751c71fcf42ff71d5e975279fe6d0d1ee051b167d5e23b86f56b4cf99b5c52a8b22769b9b7c3d7d224afb3cb3e77d89837765a983c5cddc1e3174f55fbe955a3a4b13d87d9746edce0c31d1bac44eeca2cf0973c240c5b0c";
    let A2_str = "3d4cfe9cbfe6987fb8ea7a04f62c9ed28ef9f8ced54ed453efcdaa3621d00b77b0e7a7dc12b1d582032f0b98e174aa8dcbc514b88f39988a54e6ae6be1cf4f8746a86a8837fa6f10fb999732616484422cdd7ea7e501b858f6ac3f4cfdcfcd467301ec08817a9a2955e274391d16f57d0667eadbbbe77c7b48509ca21b5ccb254d879f13b10a2dd3d089aea44ecc8bb5efac9764f86630460a5a9ea718b4394dbb091ffa5bc28624f6004c4068a7b9dce26d84d02a5a0b50b1bc168d7b9ab63b920a1ab0f7b35c09e6ed717331281699648e4f723916033dc3784491c12470bec626c3a689b95b746c342cae065dd8a8592949148d0a4a38ce96238cb6c6a5bee6e272e67462a0afa95c4f06ad7151ae1742483d1ad3203718650d40907c717eed90dd7b4e3479b623b973cc1c981f57ea2c24e76a13c982a184a4252cc3e80d0d5a5e15";
    let (A1, _) = Fq::decode(&hex::decode(A1_str).unwrap());
    let (A2, _) = Fq::decode(&hex::decode(A2_str).unwrap());

    // Kernel Data
    let P1_X_str = "ab710145f609f21b96a9829e5e9ca411431f13fe7155e1a6424013967a326a691c701d3fe90cf1536cf504fa70e0c49fcc401b7e104bec8f5e648b64bf54694c83377ac40e795be89f394e2c581a79db88f9d58d5c3114d4c07488ed95467747ac88add4ddc2f8aa95eb47be54c115b294b34fa8fb8858224b5289cccb1b1fbef87cbb7b9f1acb012c9e4341b40e7db00bdb4ffcfe88a27175aa488e887c1cf1880cdaa0a1853d95790a0902169eee111a7cf76765d92f803f615a3fe73fc01d9d245bb0325dcb5f13017c344206b583d3324fd466509308f329c3e40a6ed152e01d763eeebc2525a3e1056f82831111b54684ebcf9776d4f0be11a92de06b5c47fca1188579177bca071a152e16133e17a055da98e6bc706557f35942c0e3f665a2db4c489ac91cbe484b3092fffba9b545799fcd59626fbfe5fb77c47aef6aeb5b7616";
    let P1_Y_str = "c05a766bc0b28b4977118f646c78d43a26c1f8b7ae75847f63f6ed94e4247d581e65e85a4d86eebbe8f7928670e8962e3c51cd843b7aec237c5e1d1a9a06e7a74c7a3672e0fd8f12704bb95ec930599e5dbb2790c9480f345feb74acb27050427107aba7227730e0a6f1951a3ff22e459ef26786d93778d6f464c47ff85773985726614a91f5da610a3cf5aad29f59bd70a7c3e41861eaa0e3636df16ec987f27816753310e1ad9014a3887b2e86055cc0b087fe71148fbef4ed961df36f7a74920904709290ac3215f35056257fb284c08b02516c6ba3f7c5b66001055f0f606fad2c24b4ab975012b10ef589e1600e90f02decefb45895aa28f70d931ccce1cb17a7a5383189e10927ec2ab1dc3fb159496c7f947fcd8f6159a04ae8fd1d401929978ebb0e7d2a7ea554cfb9230351d3b39f9da0975477e17ce9291f807fc39be24813";
    let P2_X_str = "495a0140eb8101aa90f13bdab5c65e53a655ef2b8d2280b14396f498dad675305b822928d54d6f81d494d5e639e28d6c29585e788454609dc69653e5c1294b6ffa69c17e5737a8bf506c0da45a792140d167e8d67fe0b39288f0cd67cb33989f147a8d88189967e384eb54d40b105424eedd368b89b8d26beac0ebbdb2fc2bb12aa105116ecd2d076fdea288f15af4d89f543a85077db767e2ca19c7c94f3331330d5d888cd541b6297558e74cb8da409ffd041094b0e8f6725fdcda3052dae881bf1f53bfab770946cc8167ae8dbd0c5352c4153a64f2e9999bd31ae1a8ce27ad9816d5dea6bb7101a856c3804e3d26eee3b5656c58ba43a00460a9fbd4f808dd21e998abe955513a11ee68e2bcd80769ed5a1876e2f7ddac97b50d790b494abd1aa733ff76956ed56fe2457c9175a8504a52bbc9ce1ac21b9a2c2b6fc517af15b59307";
    let P2_Y_str = "eeb2858418f1768adf8e61b803381587a72dd971a982a2cf40cc263905feddf6fd4131a134c89bcc46f8a758fb1ec03dfe4d973fbac5f2533774adbe8ef0be3b7efc4650923e36f79d906d443ec7fef2e97e1dd92189ebc0e2a29556dfa821508b353b48e8b231785ba7c95e38b9ad847b08e2949be98eacfdfaa3f7837b5c436dc3034d78b6d975bee83943981c25cabe5f33a6a39faea86c4501989cd51516f6033be07b5580b03020820669a0d84ab8456837f90c665526b41ec7fe76d3a4e0c8133890845b5b92532cfdd76c1171b2b824a9bbd7b5e1895798d313f2145abd944b91d7a0d4407b608584e72ac91e4d24472f670430bc7c7bf7c6e9829a4063a254d6dd245c52fe833494b430e4cc01f34a50d823a908427908f7ddd7704bbd0ad34bb7e3ae1d10061d214d6020ad2390730393ce8fd9da5afa8252923b88463ce20e";
    let Q1_X_str = "f1327865b3ecd10fd9af8a10a7301c031637d39d66da8caa9952cf5b9e1ab8b518303a2eb8dc913bef7fb352242e3d80b27056eedded2bc1fcf5e70bb5afd610ca066c1e34c64d8801d98af75c459f805a03719a7ac29ce07b5e5625a8eabab6c705b87d217fadee068dec01c9adf4071f5eace53cc7a50914af62dff0d5343c5efb28f67fa3f4018a7192eeb0ba11bcc034f7425bce1adb07369c92268a1f3d170547805c63ef9f44f099cf3f1733304fd012a3feec1ba33c7f9626e24c83ca089ff82005e39724dae7202ea80c5e3d1271b2aeceb6046dd26d01228385f8fb4493c05d0c0a6f9c3eb9fbe4ed0cd7dd2e9bf0896e6604b4a2f9e3b4ce8b7eb8330002f497bba8985854b8652b29a55b41af08ccc4e3db74d71d65ed7ecdef8ed821ae3ac65bfa60f23b089b674d21413a6bfc04cf0345ad6ee50d16a65b84e83564df0c";
    let Q1_Y_str = "6c3b1385160d8c41e4fcd4b0f12908a0ad58772a35b34fc0a15159a5a164ebca6f53873d8b44ea79435e1d3c7afa318286e2c952070f81cb46d8e55a9ce221fc0637e5f3e5dba2422e4858e3ce6720f5070e7ff93837c2d1df96327865497dfde7160e7e6fce952ddf6a4c355d0c65127840d73999e7b81a447d07e694d873639957fb67c8d3f5e10622a03e151899ec68f24c1879fa8f9e44f4749014e27990e60272fa04158db4ad5b7274eb88a456ed34d5e67cf1eaa8867f2158f85355c5ec68bb04b086622d28ef2ff62dfde9988d79e5444edadb59fbc999c9fdaa33488c94b178d9c218109682ee4ba51b8709fc81bb66f0683e219b6c5b9f2a7f05dda514329645d7858fa65093a86fe8864bf1dcabe9f5d3c952c27ab45a4491d9291f37b2876bb62e9d96ab3fa29bb5133c18ffab247e024cd30b23f2111eb7d4a3a9759b0a";
    let Q2_X_str = "34678391557c9e2776dff7177e4c6fd8e62442ecc0aceea3a2b7969b9720c493fb4d96affe30264e5ab4c5e73f63d35da8ccde60fe738470c54b849fa9a99d78a415e651c583f18990ce9483cc905c2338dd704b738c9158febea12db9bf2c546a101a0fbf6e4923f9fd89090ef4843a2005dab3c0a14f080b7589cec37f7ee9477da4cdadf3311ae9a599734ff2de0e9fb364627eb0f5ce10ec05cadd9f8e4fa904b646f943c1cd52a1ec46145364f1edc41700525db027b0014b97caab13f747fa72552e51e41e185064178d643ce49547eff8c06bb232ab372f6d4ecfeb8205ca7c42f70f0e8085d15a91d3801bdec8aee5ea1c516410df8fd3f0f09be42876614fc25ee675fc743d1d7f58b6a31af37370a5002ec639ed014f4a2c6fbf1b646b6a24b6cec49e3c9f7e1e1095e32f880640f47ed835c50aa7ea313d9068f1d7e6980a";
    let Q2_Y_str = "7c33493b0c0ab5bf88131c9bc9291b428d6812683d1f82b0f6c41dbc95e999b336a0a9efd0050255f0f290d76291ea7c72bb03d3e3130a2101ea5b4dea3a62397f56242e5ace05bf2754db2bd51ad90f37d92a4b7fdf9e705b4788b2a9e6088faf446475bfeb4cc8fd19b6655c5953611e673d3cd354132cf484e9e395eb557991ddd4d34a63036ac95801769d2ded726c8c77ae8d9bdc29dc42c34a75d03e0489047870e6f93ebcb5ec4cdf849b2c7f28404f3380c9b31b02e43e3c8a1c5e892aaa68d1c260f0d903055021027db1ac829ea2603f775e9c3b6d5594accb05c0656d9dcbd9c030510d1c3de017cce56e16d16a4e8a0aa5f31e1d9fc565cedd7757bac076f4d06fd9b24887af3d5829a0df783f0485890f9c1c1446cb99dedd1ddb81beb6f48f8c6a1928aa053d9c3c634ea104388d28a592e4c2aab38fd01999aa797c0b";
    let (P1_X, _) = Fq::decode(&hex::decode(P1_X_str).unwrap());
    let (P1_Y, _) = Fq::decode(&hex::decode(P1_Y_str).unwrap());
    let (P2_X, _) = Fq::decode(&hex::decode(P2_X_str).unwrap());
    let (P2_Y, _) = Fq::decode(&hex::decode(P2_Y_str).unwrap());
    let (Q1_X, _) = Fq::decode(&hex::decode(Q1_X_str).unwrap());
    let (Q1_Y, _) = Fq::decode(&hex::decode(Q1_Y_str).unwrap());
    let (Q2_X, _) = Fq::decode(&hex::decode(Q2_X_str).unwrap());
    let (Q2_Y, _) = Fq::decode(&hex::decode(Q2_Y_str).unwrap());

    // Points to push through isogeny
    let L11_X_str = "dea2bb5f803a0281b2c1dc948ce739aa77ccb66d167b62b069c4e6c5e8f2d6a21bb7b385d134f04608caf3f129aa29ed6ca8be9369ed27c43535bd282878c1e8c8907581110a22df7b2ab1358b8aa50fc665b7ca399aea781ef0ac250277583cb97dc7116f98e6a45810c468a9d5da4825f3af65c4a6b96debe9fa0b5f7805312eb5403a7fa97da650e475d1daa6a93db1e5bfef4ef4e977fae52f416a2af1cf2c04810a0c60deebb39bf19d01e74e0ed02cb85518d7285d804b56340e7e68bd3e59312c2b4cad5ece56415b5ec824c4be694df9255ff04cb3905af56adf4d6473efda1922ee615d1d3c0741eadf41d1c9587d49ceccaed909b8b0174e61ce6f0df628a03da0def573fd3b7d4850bfa751b73a37839ebb484f9a9f112d95d1ee2522dc3703c5102d957fc3dfbcc24d5618a6cf5b106ef1f23590fc603699aa5887f01302";
    let L11_Y_str = "abe10e6a88932185ab71f34c0a808b1a6e34844f6bea614c04515b2755f84947d13eb29658659a4d4c66a0e2031964448f69b14a3cfe8a08533981747a82d1e648bfd9faba0b5349fe1a057e356b4fc9fa6d71078173665cfb07e19b03e2da3c45aab7e70ca1c80a0a9e0bcb2b551103c58da65277897c2470b70a05595ade7f7950c04f57497c8421d9081591529f54942fd749fb86c28ebedfd85e5d3e8f92131275c689849076add344565b166beee4d889791d9072f313305ce5c8c1f39c6c339e986d852672a1eff6a7f08ed7c99abf884db79806adf520c2fc8a34d6ddffa600454f181b563375de64f62e16ac48aa413301d5d4443ad325552d3381c4b18e317ed659978967a71d3a68c9238b6828d4f9483f5c298a872203036e5380112b1da86f0f7b44c205d20c1290eec4c644f0b5aefa1dbaac61142bc8d4232ac10cc30c";
    let (L11_X, _) = Fq::decode(&hex::decode(L11_X_str).unwrap());
    let (L11_Y, _) = Fq::decode(&hex::decode(L11_Y_str).unwrap());

    let L12_X_str = "67bd25d443ab35439cde45a4cae9deede7c7ab808d3cf080ae2d2fa8ddfbbde33fb472b24d0838a80ee4af7f4d3fe5ae85a1becad033785fb7024b4d203b3372c9a71816f0f16108fe4f9f989dec5d8427af9f49a7f3c43c3f7e8eef91bd8126cdaf3e953ae90fcb679e9d84beca4b803af9f8f6ab0af82ab00693a258121450afa5fb59227855a0e6b7301362382b8bd040d72d5e48a5627bd73d1b8c4c3e94d00fe3807cf8708a8e5752c0c4b42e9b68cf6640bda42a05efb996560aac1c407a6c0b36bbbcb9c6355cb897bf7191aafdd8ab7582becd10a5f57e40851debd9076ad86d4c425152e8a5609ee2a2c50737d05e09a2089cec0066d6e8a8f8e613aa06cd268e09725ca0ea0ad865753779e2350f7324a41145f57eec8ca3527fcd28024cab1cf831ae25398f88ad3052cd6a8955573af4e79e77eea6f3353282234072140b";
    let L12_Y_str = "fc5fb5ec8552da8cdde3fb67034abcd9832e2ae0ee73a46de6d17d1f3565eeaffb053011a34e951572ac8ba83840f713c104f171fd76979ccdb3fdba531df2fe172d37bd01a7246a3132887c6ec5cc76a6e8a031df962a86bbbdfc0ec813b501a1246a34d784a4c4bca9cb10785a7eea9a068df98907769c351292287cfa6f441de196f2406e966ab412140789260f8c918390c8149df0513bf3fda0eb1d275dcc00d30f507ce8ad17a91cdb083227b971162b04284eee36506ccf6c71b1ee07d3f2ac06e74192cd560e3d6093a4375b18b7b6281c9e1cc0583495e5169e85858ac2be3e600c6a8ba53f6c7bd9ede931fca251fa0fe2a218ae3d47565d45fdcdbaea8a262a22f9b472ee2c6fb3df9963fa2a6d8899485659ea42110f0b5dc8c92855de8dee7e64ebae6c45af01f83e36a7a3f8775808867cb165e04918525889d1661706";
    let (L12_X, _) = Fq::decode(&hex::decode(L12_X_str).unwrap());
    let (L12_Y, _) = Fq::decode(&hex::decode(L12_Y_str).unwrap());

    let L21_X_str = "5b3bcf51421852af967115e3824bbeffe9b9aab16762e530e35c32a77b0b37c62a7e0ce96f7af6b831cf7193f807ed48749ebee12c5fcf8a2fd9bf187b5a58856888f99c81fb89644138d0b136589aa28da0feb085cf3d2a43c0d54bd8c87aec9b818542ff96c0a95643402026f29eba75c9cdf0deabd841bc1f1e88f5dd04bdac64174c5b7dc7af8428f3c74c828d5dd37ff1afda9e00d4a4b8dc3d579bb1425809d8a751619517a25f9ac408ddc3391c73b170127d1e2492d8877ad069052b6cba9159ff2590bc8e628496002bca3665174d414ed67d842f93931418560264ebb0359f9ae42ca64ddd75834cdeae6e7e7172f981a419406fe24401c4d70198cfadb6f301125b7c3f37a65babab3f347c935d0cbf8ef8e7dda642767960e73380fef72976e386007bba68fefc04ab876a50db3e0c7f9b12965b61e4fbde8f887c47f70d";
    let L21_Y_str = "5ff5ef87185ee6951c538ef7f5e16b17c3c9a28b8e7828743e45472482d01062d5d616cfb34d86fb0a66cf65aa7820353b7d8df025e91b6685009ef4b2231983c3b428e79e11afe00aa7f69a7b4f4f1c09f51953b44ba84b2b0d5e80c47f66a5e883353acca8b79cd65b5c77185c028a5af8b08cdc33b01ae35be1866b932cce8ef2803800a50f07dd6d835503102ed21871f13814051ae73d8d1e29ef85a94032082996a0c8e58c5460b9dfef65311511e8030bb434ce7e4bd98ee7981143f2ac175e141776e583ce6ad2e2e06bc6ad2e20effe7969b9cd57ac1adb7cd89decce7347fcf1522dfed1ff2663d3942fe9fd86dd798d9ff6edb0db17472259e5b63086e4e8e04f112c474f8a4b64dd5a89cfb176ea9742f954a507b21dc3163bc1609dd5077daec60c14851496d3e56612210a8bff825890752462c02e8b720d909b685308";
    let (L21_X, _) = Fq::decode(&hex::decode(L21_X_str).unwrap());
    let (L21_Y, _) = Fq::decode(&hex::decode(L21_Y_str).unwrap());

    let L22_X_str = "4eb3160cce9bb8a02fbd87d6f070cfdab2f0465f2ed692457f56a9699f3d2a51385d1b0ea6c112676da50db2cd916295aeae5a0c8fab327c9019ecdbbae67d0677c38688b7988a99bd07bc4f67f47de987f535bdf4f7af89b6f17aebd5841eebc675a3bd0622b4445a86f469ab6478760e8d54f4810e60f82f25df97957768708b2373e28c5bdc051f5eb829ff63e25f024f546d7c0d2d30a90d235961ff861f14163fe43d8b8c9b8fc5f8bc8ee7970f3aa23ab78eb4a9a836ff087806843a4935860aab038b0809f7df120f5eeaa5fcf6c936158448cac458ac9c7ca3e89d830efde790240a79a515c4cb76e6e93ed2b121bbb75397153154ed16345358a9b2e6e28c70cef892dc01fb13ee017d04ee3848cf67fc8107b80b02b139ac9cb987c47cc3b8746a49d9e21e68643861d94cfbdce8aaed2ec3fd81be6b4aeefd27a3d85d7c12";
    let L22_Y_str = "6eb74be4d694cfd39f3a216ea77514079c66b8dff598af6588683946493a691d812685e352c218f95f465e2b21da95012c409cbf5c9925ae5e2c3ba6249164641d0314f03de05f3f9cd6851ae2681f04a7d4d04a76999ff6ce0dcb73f2d1e76fe9c86a403c62f4de281c3bccbfa72410a9bf26a7d2959577be8d9b23de9c634fd096b468e1da9ee98a4a449ac99bd590c1f9b60615dc748b6d7966e61e373edb7b07c97cffcd09ebc8daa380977a10555f800229e2bd158acc311e4e0bda66dae5f2b6d0745d62ed3d9a7d5415ba6aecfcc7db6eeff54476179c59e5cdb6f62c32f634b2de573181629b25d3febea4943982f72b535389183118caf16685fdb5ccfa73e609f9a7a98b7b1209231d0d9004cdff6b425cf578a4f586473e9a4fbf32b6fbae0cd1ecc1ea48d6ec950b2c42fe9a841a9cd559ab0154dc60b73ea290d1d37613";
    let (L22_X, _) = Fq::decode(&hex::decode(L22_X_str).unwrap());
    let (L22_Y, _) = Fq::decode(&hex::decode(L22_Y_str).unwrap());

    // Result to check against (Debugging)
    let im_L11_X_str = "5d1e75b96024eb526a1a8ab9ac93c105b12107125965d5005557a5255077a56e4b435a66fc17d554cc6b0639811ede4e7a06f4f93bf2ff9858115305eb3895913a737109c7b3e446cc66a7fe4b6a0836cb46c8c39d542a6fcb7d7585c5753196b80f2e53c8f385ff95dbf841303818fd95545a9e7bf6216cf2fac12d544827310b65c9c729bfea6c8ae27b8ab5947e183656eb8391f0bff1b983c1f9f317227d17046a442f190320cf651d95b3fc3537dd5b21616e3218cd8706f861a35a8f09dab19205c574b4b3ecfd099238189ad3ab900928dd43641c4310a45715cc836d76d6863e54847af02a5dd9c29adfb067b0c7f16f8a2d2d35c64427b8606e62e59267938eb9e8d1047b3676aacae48d55af42f0fd5a250828b88333c7a4c87846b9a22d489d6deb3746ee7ae9d173cbb4bb66071fca218cbab6c64b3e2b6fe230efcb4e08";
    let im_L11_Y_str = "c25bd83bb502e6a8ffdbafadd577c2529b457cdf1d470b764b61f05a0a5a331c20642902f1709d48295d691e8f6b47b423e653ac2054b59c04f980aa07586b5e24213dd335e7262b91c4c1b620817802969238a1f459b7f9bdc2ce31e5a38d5fce542132dc234c1645d95f171e7aea4a216a11b9eea3ee655c4e35f71f74a77ffa84229b215f4e51b8381637cf04817225f6b2333aa5de78f8af4f8d757883847e04d77ece30e000a72050c4664fd01a7c8e0c8afe4a4c2cb3ee56fe187effc9401a275f88969f7649d30bfcd1f23b64521b7dce274cca02289c3e41935d87a0a3a52d36342d617f2ae63cfdc0cc6ae2c20f11a61da9b3c695609a2ac4f818e4d78766597a0a44fc8b8428eedd62839a434318f5edc46db5279f2fe5794b89523ea97b3b7b532ad5b5c93d948ac175ab2131996ebce4d591408daa2296fda1210229df01";
    let (im_L11_X, _) = Fq::decode(&hex::decode(im_L11_X_str).unwrap());
    let (im_L11_Y, _) = Fq::decode(&hex::decode(im_L11_Y_str).unwrap());

    let im_L12_X_str = "b2b1c1489ad8010659ed45acee631ecd4f1654dcd014822d91963f53829e03dab4016ca088d377d9782cb983e78888b3f6b2fafbf9895d114dfc115899b08a4accc4ed467e3a80e2016542764f03afc3c0425c729eea94cf15738d1c50b68bcfd753e6f19a5bbbaf63abbec2b4ce10d29a255f9374c8bf718cf122e5b336c35872045afa5ca7c5d00ba1fc7b288bb3e8ff7729c2a66c97612ac497068324c662870240182b54358a2f3144c770800d05a4184ea4760161dc75a6a7d7604536dff6829af7578e8eab2eaf439eff87c69fd9d10a0ed06f34aff46cf15de71e6beaecaf7449b2e68a9f94f49b2e4fe98cc12c37e445b02762ab4fd3e8f75ad4ba1026da08d16f989b2a9cfb6949961e45b389ecff7406fb259821c2cee000fb1bde62619ceab4d3c1f4943c03f24a22ea0123fab70118c8c2da9d7540848cf7137bb06aef15";
    let im_L12_Y_str = "aa9e26c4724b97d8ce6a65ba2dc53d635c0105c6d1a7f6be7f3838714b5cd52b8c8c6d5a492331fabc513e1f759dd7417dcf4fa44b3354963ba32c6a2e3c22d0bd03484d96a980a9390acbbfb3243577f628d37865dac091b13b8c8f42b264a4b0ad3c1a5b6ca892f13eb3750aacf234e8db035ff6cda25b524a62a47b84da3e734f72e8aff2370e60ed2b2749c6d77085ee83e843bcd018660ae5e7458f1fd98b021559ea8727e553a4fef2200a9cc30c70d9a5aee5d8c7c35e0f9107ff9cf964e8c2d38fe56f7db19d067fe70fe4bcc50ad32e7034c62094731bffff7207e73362cd7be2048cc31242052f9fc77898d48ff8e1fb289720727d1594e4897edbabf8dd10c3d1894cc779483c12b3633341b1d464c58a22d21dfbfd8f24ddd7614c95a8d99ec8b6fcb0442f51fc3930e8944d15da2a8f588c36f61bc8c3da4a2815d17e11";
    let (im_L12_X, _) = Fq::decode(&hex::decode(im_L12_X_str).unwrap());
    let (im_L12_Y, _) = Fq::decode(&hex::decode(im_L12_Y_str).unwrap());

    let im_L21_X_str = "7b0af2b62ad1418d0e6f45e4f5e1e23aba8b21e603c1d327ecf65292fd1361210388c532d4076e3134a94ca28b936726f10a34062018ecdd278af21dc5aee3fe5129ed7592228bfb64d17ba67bb3dc02deefe95ce32e143c6e899eec2bd5b860364c9ca83fb06290d7a1fced8070049da7a2f61a279d79a2a6f3e3fc86efca46a8291e6b263c4df0cc086fa3db09584f321db9ae9c63d343657d3c8e9e2cb3186f0985a583085a80ff7c987d49e547c1ce65e99384a62a5b600939dbecff800962f22a508bccd7c87ce8fc0184fc75106534cd5a55f352a57b3cd3e55141deb2fc141768f2bcb667178496a0d89fd77947cc0d3fa0541b109fad7f39387f3eeca9ae1cd8cbb70efe29f6913dedd30b2a136f3660be6ffe1de52ea9da31ac9d392b667c6ff9feb5bee4c1b6da31ea62b9d484c7e63593324fa98982e6f2b4940a79471210";
    let im_L21_Y_str = "34652801a908dd0b5567b615b787a45b865b6a8fa4a5aecfbc7b434e671204705ee2a84de0c8ec82b60cfa259b6cbdde8f76a777e0e347bdd7c39c0cde925e61834fb6982de4d0cb1fa5b1144989013ec274b4b59878a0475a6dcf7d82d3b3104e2e2e7a01d5d1d308b705601c38368d3ef2b447c9eaef6629283a34fe157fff1e0392b3553f4da2a39ab66401d197239d897ed81f4d327ff448436ed2383a72f7121fb6e2674ed32f0b245dd96f50b5090c7acbffd9132da9675fd09322d74db9ffb487d5a9d9a259519cbd5d4f8d5f6553258e2f66ba1e8aa7bf543c8ec1ee32a2038d861c26b908a2de5ade9c60de6cb91cdddba7eda3928b129dcebb1d1976e3fc1a505baaf0bee43bfe3ca6bd826bfe37166e1d777a7310b8c58e7e136d894e730b77a68db35d3e026ded76d7bd80f9ebccc89ff0361084a8a2fc74c6dd01ba2c15";
    let (im_L21_X, _) = Fq::decode(&hex::decode(im_L21_X_str).unwrap());
    let (im_L21_Y, _) = Fq::decode(&hex::decode(im_L21_Y_str).unwrap());

    let im_L22_X_str = "6f8c129e738d35e908803a47aac5f45e0b6bca7ad4b0c59325d5169a1151ef8273f38e90e648cd8de58487544b8f2b9ae73697ba60d7ba9220bb25feffebc4638f7a7a8ded1f25367568e96d8c1e2bb38d4a9a13bd0ab02785695696a7d63a11f89622d2b373cf8c2412129812b30a62d032ffbcc192ed9398e1616e61c05a186d6ee3b3280ffb839387e7f4e0c0180a468cef504a39709e7ca1ed0172d0dd3655120857fe033936490b7c3ab0181a8c6b718c9ba86d36002c18dada113c735e5f628a4c0644bea3f78b16c2a9c16a5cfcdd5859827b5ce2c2752915d6c73b58f709ea3e52b695adca6a9c36ca9e4dc655145c9cdd2357e2cf969a5ad8486869d8e9357a0465345c1992a0117414f1dfc7138c5538ee7f6da349a8cbca8e599dc24064b5c162220bea04582fbe528a1f0abe9bca9698305c1f63054239f4dbfeb490b002";
    let im_L22_Y_str = "dd0683896fcf1abb06040dd596d5c25cc993a1776cb2d99f56c7212f029b8f68feb90cb54b3340396c189d42e619c2af1f42cee9a3c6ca7d5a35e7958e45bd203c6d6e96438c0a887e2ac0ad50e17a0ee24bc028d908a7edfcb9e9b9ff4d03556f2163090740ea97eac5de8d568f64ec6e7fac381b148b9f7570cda0ac71b01b0b561db179ce3a5adc29b9f4d0ac40d110fa01d7586f852296d99fb91d6feeea110c2d1a8499b556a372fbf5a7dee4d49e264fe452ba0ce592196957b9aeb226e7d053c645bc2a5c86b4955fb9b73bdaf480556532a58ed756f28630124725f8ef1b7fbc80f835175ba48baa9d39ffe9fc6d534aaa3d1631834f13a78cc1b356850ee4be86609e8eccc6e0355aa7ea26ec1f65cb5e6bf7b84a6ca5d76e03e72e59744bb1137a2c7b4688f0eeffb28360430f32589aa3fafba83a43cbfccd31a9645ed610";
    let (im_L22_X, _) = Fq::decode(&hex::decode(im_L22_X_str).unwrap());
    let (im_L22_Y, _) = Fq::decode(&hex::decode(im_L22_Y_str).unwrap());

    // Curves which define elliptic product
    let E1 = Curve::new(&A1);
    let E2 = Curve::new(&A2);
    let E1E2 = EllipticProduct::new(&E1, &E2);

    // Kernel Points on E1 x E2
    let P1 = Point::new_xy(&P1_X, &P1_Y);
    let P2 = Point::new_xy(&P2_X, &P2_Y);
    let Q1 = Point::new_xy(&Q1_X, &Q1_Y);
    let Q2 = Point::new_xy(&Q2_X, &Q2_Y);
    let P1P2 = CouplePoint::new(&P1, &P2);
    let Q1Q2 = CouplePoint::new(&Q1, &Q2);

    // Points to push through isogeny
    let L11 = Point::new_xy(&L11_X, &L11_Y);
    let L12 = Point::new_xy(&L12_X, &L12_Y);
    let L21 = Point::new_xy(&L21_X, &L21_Y);
    let L22 = Point::new_xy(&L22_X, &L22_Y);
    let L1 = CouplePoint::new(&L11, &L12);
    let L2 = CouplePoint::new(&L21, &L22);

    let image_points = [L1, L2];

    // Points to compare against
    let im_L11 = Point::new_xy(&im_L11_X, &im_L11_Y);
    let im_L12 = Point::new_xy(&im_L12_X, &im_L12_Y);
    let im_L21 = Point::new_xy(&im_L21_X, &im_L21_Y);
    let im_L22 = Point::new_xy(&im_L22_X, &im_L22_Y);

    // Length of isogeny chain
    let n = 632;

    // Strategy computed from strategy.py
    let strategy: [usize; 631] = [
        233, 144, 89, 55, 42, 34, 21, 13, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1,
        1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1,
        1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 13, 8, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 3, 2, 1, 1, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8,
        5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1,
        1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3,
        2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 55, 34, 21, 13, 8, 5, 3, 2, 1,
        1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3,
        2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1,
        1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1,
        1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1,
        1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 89, 55, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2,
        1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1,
        1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1,
        1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1,
        3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1,
        1, 1, 3, 2, 1, 1, 1, 1, 1, 34, 21, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1,
        1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2,
        1, 1, 1, 1, 1, 13, 8, 5, 3, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 5, 3, 2, 1,
        1, 1, 1, 1, 2, 1, 1, 1,
    ];

    // Compute chain
    let (E3E4, images) = product_isogeny(&E1E2, &P1P2, &Q1Q2, &image_points, n, &strategy);

    let (E3, E4) = E3E4.curves();
    println!("E3: {}", E3);
    println!("E4: {}", E4);

    let (P1, P2) = images[0].points();
    let (P3, P4) = images[1].points();

    let P1_check: bool = P1.equals(&im_L11) | P1.equals(&(-im_L11)) == u32::MAX;
    let P2_check: bool = P2.equals(&im_L12) | P2.equals(&(-im_L12)) == u32::MAX;
    let P3_check: bool = P3.equals(&im_L21) | P3.equals(&(-im_L21)) == u32::MAX;
    let P4_check: bool = P4.equals(&im_L22) | P4.equals(&(-im_L22)) == u32::MAX;

    println!("First image matches what is expected: {}", P1_check);
    println!("Second image matches what is expected: {}", P2_check);
    println!("Third image matches what is expected: {}", P3_check);
    println!("Fourth image matches what is expected: {}", P4_check);

    println!();
}

fn main() {
    sqi_example();
    three_lambda_example();
    festa_example();
}
