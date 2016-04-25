
function [M1 M2] = CalcP4PfCoefs(glab, glac, glad, glbc, glbd, glcd, a1, a2, b1, b2, c1, c2, d1, d2)

	% coeficients 
	cc1 = 1;
	cc2 = -1/2/glad*glac+1/2/glad*glbc-1/2/glad*glab;
	cc3 = -1;
	cc4 = -1;
	cc5 = b1*c1+b2*c2;
	cc6 = 1/glad*glac+1/glad*glab-1/glad*glbc;
	cc7 = -1/2/glad*d1^2*glac-1/2/glad*d1^2*glab-1/2/glad*d2^2*glac+1/2/glad*d2^2*glbc+1/2/glad*d1^2*glbc-1/2/glad*d2^2*glab;
	cc8 = -1/2/glad*glab-1/2/glad*glac+1/2/glad*glbc+1;
	cc9 = -a2*b2-a1*b1;
	cc10 = -a1*c1-a2*c2;
	cc11 = -a1/glad*d1*glbc+a2/glad*d2*glac-a2/glad*d2*glbc+a1/glad*d1*glac+a2/glad*d2*glab+a1/glad*d1*glab;
	cc12 = -1/2*a2^2/glad*glac-1/2*a1^2/glad*glac+1/2*a1^2/glad*glbc+1/2*a2^2/glad*glbc-1/2*a1^2/glad*glab-1/2*a2^2/glad*glab+a2^2+a1^2;
	cc13 = 1;
	cc14 = -1/glad*glac;
	cc15 = -2;
	cc16 = c2^2+c1^2;
	cc17 = 2/glad*glac;
	cc18 = -1/glad*d1^2*glac-1/glad*d2^2*glac;
	cc19 = -1/glad*glac+1;
	cc20 = -2*a2*c2-2*a1*c1;
	cc21 = 2*a1/glad*d1*glac+2*a2/glad*d2*glac;
	cc22 = a1^2-a2^2/glad*glac-a1^2/glad*glac+a2^2;
	cc23 = 1;
	cc24 = 1/2/glad*glbd-1/2/glad*glab-1/2;
	cc25 = -1;
	cc26 = -1/glad*glbd+1/glad*glab;
	cc27 = d2*b2+d1*b1;
	cc28 = -1/2*d1^2+1/2/glad*d2^2*glbd-1/2*d2^2-1/2/glad*d1^2*glab-1/2/glad*d2^2*glab+1/2/glad*d1^2*glbd;
	cc29 = 1/2-1/2/glad*glab+1/2/glad*glbd;
	cc30 = -a2*b2-a1*b1;
	cc31 = -a2/glad*d2*glbd+a2/glad*d2*glab-a1/glad*d1*glbd+a1/glad*d1*glab;
	cc32 = 1/2*a2^2/glad*glbd+1/2*a1^2/glad*glbd-1/2*a1^2/glad*glab-1/2*a2^2/glad*glab+1/2*a2^2+1/2*a1^2;
	cc33 = 1;
	cc34 = -1/2/glad*glac-1/2+1/2/glad*glcd;
	cc35 = -1;
	cc36 = -1/glad*glcd+1/glad*glac;
	cc37 = d1*c1+d2*c2;
	cc38 = 1/2/glad*d2^2*glcd-1/2*d1^2-1/2*d2^2-1/2/glad*d2^2*glac-1/2/glad*d1^2*glac+1/2/glad*d1^2*glcd;
	cc39 = 1/2/glad*glcd-1/2/glad*glac+1/2;
	cc40 = -a1*c1-a2*c2;
	cc41 = -a1/glad*d1*glcd+a2/glad*d2*glac+a1/glad*d1*glac-a2/glad*d2*glcd;
	cc42 = 1/2*a1^2/glad*glcd+1/2*a2^2/glad*glcd-1/2*a2^2/glad*glac-1/2*a1^2/glad*glac+1/2*a2^2+1/2*a1^2;
	cc43 = 1;
	cc44 = -1/glad*glab;
	cc45 = -2;
	cc46 = b2^2+b1^2;
	cc47 = 2/glad*glab;
	cc48 = -1/glad*d1^2*glab-1/glad*d2^2*glab;
	cc49 = -1/glad*glab+1;
	cc50 = -2*a2*b2-2*a1*b1;
	cc51 = 2*a1/glad*d1*glab+2*a2/glad*d2*glab;
	cc52 = a1^2-a2^2/glad*glab-a1^2/glad*glab+a2^2;
	cc53 = 1;
	cc54 = 1/2/glad*glbc-1/2/glad*glac-1/2/glad*glab;
	cc55 = -1;
	cc56 = -1;
	cc57 = b2*c2+b1*c1;
	cc58 = 1/glad*glab+1/glad*glac-1/glad*glbc;
	cc59 = -1/2/glad*d1^2*glab-1/2/glad*d2^2*glac-1/2/glad*d2^2*glab+1/2/glad*d2^2*glbc+1/2/glad*d1^2*glbc-1/2/glad*d1^2*glac;
	cc60 = -1/2/glad*glab+1/2/glad*glbc+1-1/2/glad*glac;
	cc61 = -a2*b2-a1*b1;
	cc62 = -a1*c1-a2*c2;
	cc63 = a2/glad*d2*glab+a1/glad*d1*glab-a2/glad*d2*glbc-a1/glad*d1*glbc+a2/glad*d2*glac+a1/glad*d1*glac;
	cc64 = a2^2+a1^2-1/2*a1^2/glad*glac-1/2*a1^2/glad*glab-1/2*a2^2/glad*glab+1/2*a1^2/glad*glbc+1/2/glad*glbc*a2^2-1/2*a2^2/glad*glac;
	cc65 = 0;
	cc66 = 1;
	cc67 = -1/glad*glac;
	cc68 = -2;
	cc69 = c1^2+c2^2;
	cc70 = 2/glad*glac;
	cc71 = -1/glad*d1^2*glac-1/glad*d2^2*glac;
	cc72 = 1-1/glad*glac;
	cc73 = -2*a2*c2-2*a1*c1;
	cc74 = 2*a1/glad*d1*glac+2*a2/glad*d2*glac;
	cc75 = -a1^2/glad*glac+a1^2-a2^2/glad*glac+a2^2;
	cc76 = 0;
	cc77 = 1;
	cc78 = -1/2+1/2/glad*glbd-1/2/glad*glab;
	cc79 = -1;
	cc80 = -1/glad*glbd+1/glad*glab;
	cc81 = d2*b2+d1*b1;
	cc82 = -1/2*d1^2-1/2*d2^2-1/2/glad*d2^2*glab+1/2/glad*d2^2*glbd-1/2/glad*d1^2*glab+1/2/glad*d1^2*glbd;
	cc83 = -1/2/glad*glab+1/2+1/2/glad*glbd;
	cc84 = -a2*b2-a1*b1;
	cc85 = -a2/glad*d2*glbd+a2/glad*d2*glab+a1/glad*d1*glab-a1/glad*d1*glbd;
	cc86 = 1/2*a1^2/glad*glbd+1/2*a1^2-1/2*a1^2/glad*glab+1/2*a2^2-1/2*a2^2/glad*glab+1/2/glad*glbd*a2^2;
	cc87 = 0;
	cc88 = 1;
	cc89 = -1/2+1/2/glad*glcd-1/2/glad*glac;
	cc90 = -1;
	cc91 = 1/glad*glac-1/glad*glcd;
	cc92 = d2*c2+d1*c1;
	cc93 = 1/2/glad*d1^2*glcd-1/2/glad*d2^2*glac-1/2/glad*d1^2*glac+1/2/glad*d2^2*glcd-1/2*d2^2-1/2*d1^2;
	cc94 = 1/2/glad*glcd-1/2/glad*glac+1/2;
	cc95 = -a2*c2-a1*c1;
	cc96 = a1/glad*d1*glac-a2/glad*d2*glcd+a2/glad*d2*glac-a1/glad*d1*glcd;
	cc97 = 1/2*a1^2+1/2*a2^2-1/2*a1^2/glad*glac-1/2*a2^2/glad*glac+1/2*a1^2/glad*glcd+1/2*a2^2/glad*glcd;
	cc98 = 0;
	cc99 = 1;
	cc100 = 1/2/glad*glbc-1/2/glad*glac-1/2/glad*glab;
	cc101 = -1;
	cc102 = -1;
	cc103 = b1*c1+b2*c2;
	cc104 = 1/glad*glab+1/glad*glac-1/glad*glbc;
	cc105 = 1/2/glad*d2^2*glbc+1/2/glad*d1^2*glbc-1/2/glad*d1^2*glac-1/2/glad*d2^2*glac-1/2/glad*d2^2*glab-1/2/glad*d1^2*glab;
	cc106 = -1/2/glad*glab+1/2/glad*glbc-1/2/glad*glac+1;
	cc107 = -a1*b1-a2*b2;
	cc108 = -a2*c2-a1*c1;
	cc109 = a2/glad*d2*glab+a1/glad*d1*glac-a1/glad*d1*glbc+a1/glad*d1*glab+a2/glad*d2*glac-a2/glad*d2*glbc;
	cc110 = a2^2+a1^2-1/2*a1^2/glad*glac-1/2*glab/glad*a2^2-1/2/glad*a1^2*glab+1/2*a1^2/glad*glbc+1/2/glad*glbc*a2^2-1/2*a2^2/glad*glac;
	cc111 = 0;
	cc112 = 1;
	cc113 = 1/2*glbd/glad-1/2/glad*glab-1/2;
	cc114 = -1;
	cc115 = -glbd/glad+1/glad*glab;
	cc116 = d1*b1+d2*b2;
	cc117 = 1/2*glbd/glad*d2^2-1/2/glad*d1^2*glab-1/2/glad*d2^2*glab-1/2*d1^2+1/2*glbd/glad*d1^2-1/2*d2^2;
	cc118 = 1/2-1/2/glad*glab+1/2*glbd/glad;
	cc119 = -a2*b2-a1*b1;
	cc120 = a2/glad*d2*glab+a1/glad*d1*glab-a1*glbd/glad*d1-a2*glbd/glad*d2;
	cc121 = 1/2*a1^2/glad*glbd-1/2/glad*a1^2*glab+1/2/glad*glbd*a2^2+1/2*a2^2-1/2*glab/glad*a2^2+1/2*a1^2;
	cc122 = 0;
	cc123 = 1;
	cc124 = 1/2/glad*glcd-1/2/glad*glac-1/2;
	cc125 = -1;
	cc126 = 1/glad*glac-1/glad*glcd;
	cc127 = d1*c1+d2*c2;
	cc128 = 1/2/glad*d1^2*glcd-1/2*d2^2-1/2/glad*d2^2*glac-1/2/glad*d1^2*glac+1/2/glad*d2^2*glcd-1/2*d1^2;
	cc129 = 1/2+1/2/glad*glcd-1/2/glad*glac;
	cc130 = -a1*c1-a2*c2;
	cc131 = a1/glad*d1*glac-a1/glad*d1*glcd-a2/glad*d2*glcd+a2/glad*d2*glac;
	cc132 = -1/2*a1^2/glad*glac+1/2*a2^2/glad*glcd-1/2*a2^2/glad*glac+1/2*a1^2/glad*glcd+1/2*a1^2+1/2*a2^2;
	cc133 = 0;
	cc134 = 1;
	cc135 = -1/glad*glab;
	cc136 = -2;
	cc137 = b2^2+b1^2;
	cc138 = 2/glad*glab;
	cc139 = -1/glad*d1^2*glab-1/glad*d2^2*glab;
	cc140 = 1-1/glad*glab;
	cc141 = -2*a2*b2-2*a1*b1;
	cc142 = 2*a2/glad*d2*glab+2*a1/glad*d1*glab;
	cc143 = a2^2-glab/glad*a2^2-1/glad*a1^2*glab+a1^2;
	cc144 = 0;
	cc145 = 1;
	cc146 = -1/2/glad*glab-1/2/glad*glac+1/2/glad*glbc;
	cc147 = -1;
	cc148 = -1;
	cc149 = b2*c2+b1*c1;
	cc150 = 1/glad*glab+1/glad*glac-1/glad*glbc;
	cc151 = -1/2/glad*d2^2*glac-1/2/glad*d2^2*glab-1/2/glad*d1^2*glab+1/2/glad*d1^2*glbc-1/2/glad*d1^2*glac+1/2/glad*d2^2*glbc;
	cc152 = 1/2/glad*glbc-1/2/glad*glac-1/2/glad*glab+1;
	cc153 = -a1*b1-a2*b2;
	cc154 = -a1*c1-a2*c2;
	cc155 = a1/glad*d1*glac+a2/glad*d2*glab+a1/glad*d1*glab+a2/glad*d2*glac-a1/glad*d1*glbc-a2/glad*d2*glbc;
	cc156 = -1/2*glab/glad*a2^2-1/2*a2^2/glad*glac+1/2*a1^2/glad*glbc-1/2/glad*a1^2*glab-1/2*a1^2/glad*glac+1/2/glad*glbc*a2^2+a2^2+a1^2;
	cc157 = 0;
	cc158 = 1;
	cc159 = -1/glad*glac;
	cc160 = -2;
	cc161 = c2^2+c1^2;
	cc162 = 2/glad*glac;
	cc163 = -1/glad*d1^2*glac-1/glad*d2^2*glac;
	cc164 = -1/glad*glac+1;
	cc165 = -2*a2*c2-2*a1*c1;
	cc166 = 2*a2/glad*d2*glac+2*a1/glad*d1*glac;
	cc167 = -a1^2/glad*glac+a1^2-a2^2/glad*glac+a2^2;
	cc168 = 0;
	cc169 = 1;
	cc170 = 1/2/glad*glbd-1/2-1/2/glad*glab;
	cc171 = -1;
	cc172 = 1/glad*glab-1/glad*glbd;
	cc173 = d1*b1+d2*b2;
	cc174 = -1/2*d2^2-1/2/glad*d2^2*glab-1/2/glad*d1^2*glab+1/2/glad*d1^2*glbd+1/2/glad*d2^2*glbd-1/2*d1^2;
	cc175 = 1/2/glad*glbd-1/2/glad*glab+1/2;
	cc176 = -a1*b1-a2*b2;
	cc177 = -a1/glad*d1*glbd-a2/glad*d2*glbd+a2/glad*d2*glab+a1/glad*d1*glab;
	cc178 = -1/2/glad*a1^2*glab+1/2*a1^2/glad*glbd+1/2/glad*glbd*a2^2-1/2*glab/glad*a2^2+1/2*a2^2+1/2*a1^2;
	cc179 = 0;
	cc180 = 1;
	cc181 = 1/2/glad*glcd-1/2/glad*glac-1/2;
	cc182 = -1;
	cc183 = 1/glad*glac-1/glad*glcd;
	cc184 = d1*c1+d2*c2;
	cc185 = 1/2/glad*d2^2*glcd+1/2/glad*d1^2*glcd-1/2*d2^2-1/2*d1^2-1/2/glad*d2^2*glac-1/2/glad*d1^2*glac;
	cc186 = 1/2/glad*glcd+1/2-1/2/glad*glac;
	cc187 = -a2*c2-a1*c1;
	cc188 = a2/glad*d2*glac-a2/glad*d2*glcd-a1/glad*d1*glcd+a1/glad*d1*glac;
	cc189 = 1/2*a2^2-1/2*a1^2/glad*glac-1/2*a2^2/glad*glac+1/2*a1^2+1/2*a1^2/glad*glcd+1/2*a2^2/glad*glcd;
	cc190 = 0;
	cc191 = 1;
	cc192 = -1/glad*glab;
	cc193 = -2;
	cc194 = b2^2+b1^2;
	cc195 = 2/glad*glab;
	cc196 = -1/glad*d2^2*glab-1/glad*d1^2*glab;
	cc197 = 1-1/glad*glab;
	cc198 = -2*a2*b2-2*a1*b1;
	cc199 = 2*a2/glad*d2*glab+2*a1/glad*d1*glab;
	cc200 = -1/glad*a1^2*glab+a2^2+a1^2-glab/glad*a2^2;
	cc201 = 0;

	% coefficients matrices 
	M1=[...
		0 0 0 0 0 0 0 0 0 cc1 0 0 0 cc2 cc3 cc4 cc6 cc8 
		0 0 0 0 0 0 0 0 0 0 cc13 0 0 cc14 0 cc15 cc17 cc19 
		0 0 0 0 0 0 0 0 0 0 0 cc23 0 cc24 cc25 0 cc26 cc29 
		0 0 0 0 0 0 0 0 0 0 0 0 cc33 cc34 0 cc35 cc36 cc39 
		0 0 0 0 0 0 0 0 cc43 0 0 0 0 cc44 cc45 0 cc47 cc49 
		cc53 0 0 0 0 cc54 0 0 cc55 cc56 0 cc58 0 0 cc60 0 0 0 
		0 cc66 0 0 0 cc67 0 0 0 cc68 0 cc70 0 0 cc72 0 0 0 
		0 0 cc77 0 0 cc78 0 0 cc79 0 0 cc80 0 0 cc83 0 0 0 
		0 0 0 cc88 0 cc89 0 0 0 cc90 0 cc91 0 0 cc94 0 0 0 
		0 cc99 0 0 0 0 cc100 0 0 cc101 cc102 0 cc104 0 0 cc106 0 0 
		0 0 0 cc112 0 0 cc113 0 0 cc114 0 0 cc115 0 0 cc118 0 0 
		0 0 0 0 cc123 0 cc124 0 0 0 cc125 0 cc126 0 0 cc129 0 0 
		cc134 0 0 0 0 0 cc135 0 0 cc136 0 0 cc138 0 0 cc140 0 0 
		0 0 0 cc145 0 0 0 cc146 0 0 0 cc147 cc148 cc150 0 0 cc152 0 
		0 0 0 0 cc158 0 0 cc159 0 0 0 0 cc160 cc162 0 0 cc164 0 
		0 0 0 0 0 cc169 0 cc170 0 0 0 cc171 0 cc172 0 0 cc175 0 
		0 0 0 0 0 0 cc180 cc181 0 0 0 0 cc182 cc183 0 0 cc186 0 
		0 0 cc191 0 0 0 0 cc192 0 0 0 cc193 0 cc195 0 0 cc197 0 
		];
	M2=[...
		0 0 0 0 0 0 0 0 0 cc5 0 0 0 cc7 cc9 cc10 cc11 cc12 
		0 0 0 0 0 0 0 0 0 0 cc16 0 0 cc18 0 cc20 cc21 cc22 
		0 0 0 0 0 0 0 0 0 0 0 cc27 0 cc28 cc30 0 cc31 cc32 
		0 0 0 0 0 0 0 0 0 0 0 0 cc37 cc38 0 cc40 cc41 cc42 
		0 0 0 0 0 0 0 0 cc46 0 0 0 0 cc48 cc50 0 cc51 cc52 
		cc57 0 0 0 0 cc59 0 0 cc61 cc62 0 cc63 0 0 cc64 0 0 cc65 
		0 cc69 0 0 0 cc71 0 0 0 cc73 0 cc74 0 0 cc75 0 0 cc76 
		0 0 cc81 0 0 cc82 0 0 cc84 0 0 cc85 0 0 cc86 0 0 cc87 
		0 0 0 cc92 0 cc93 0 0 0 cc95 0 cc96 0 0 cc97 0 0 cc98 
		0 cc103 0 0 0 0 cc105 0 0 cc107 cc108 0 cc109 0 0 cc110 0 cc111 
		0 0 0 cc116 0 0 cc117 0 0 cc119 0 0 cc120 0 0 cc121 0 cc122 
		0 0 0 0 cc127 0 cc128 0 0 0 cc130 0 cc131 0 0 cc132 0 cc133 
		cc137 0 0 0 0 0 cc139 0 0 cc141 0 0 cc142 0 0 cc143 0 cc144 
		0 0 0 cc149 0 0 0 cc151 0 0 0 cc153 cc154 cc155 0 0 cc156 cc157 
		0 0 0 0 cc161 0 0 cc163 0 0 0 0 cc165 cc166 0 0 cc167 cc168 
		0 0 0 0 0 cc173 0 cc174 0 0 0 cc176 0 cc177 0 0 cc178 cc179 
		0 0 0 0 0 0 cc184 cc185 0 0 0 0 cc187 cc188 0 0 cc189 cc190 
		0 0 cc194 0 0 0 0 cc196 0 0 0 cc198 0 cc199 0 0 cc200 cc201 
		];

end